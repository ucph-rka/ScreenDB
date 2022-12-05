from pathlib import Path
from struct import unpack, iter_unpack
import lxml.etree as ET
import numpy as np
import pandas as pd

def convert_str_id(str_id):
    if "-" in str_id: 
        return str_id.replace("-","").upper()
    inserts = [8, 12, 16, 20]
    s = ""
    for i, c in enumerate(str_id):
        if i in inserts:
            s += "-"
        s += c
    return s.lower()

def local_name(tag):
    return ET.QName(tag).localname

def get_vectors(stream):
    vlst = []
    orig_loc = stream.get("location", "").upper()
    for x in stream:
        if "stream" in x.tag:
            if x.get("location","").upper() == orig_loc:
                vlst.extend(get_vectors(x))
        elif len(x) == 0:
            vlst.append(x)
            
    return vlst

struct_format = dict(
        byte = "B",
        sbyte = "b",
        short = "h",
        int = "i",
        uint = "I",
        long = "l",
        float = "f",
        double = "d",
        boolean = "?",
        )

def get_fmt(stream):
    lst = []
    for v in get_vectors(stream):
        header = v.get("label", v.tag.split("}")[-1])
        dt = v.get("dataType")        
        lst.append((header, dt))
    return pd.DataFrame(lst, columns=["header", "fmt"])

def get_headers_fmt(stream):
    df = get_fmt(stream)
    headers = df["header"].values.tolist()
    fmt = ''.join(df["fmt"].apply(struct_format.get))
    return headers, "<"+fmt

name_expr = ".//*[local-name() = $name]"    
name_expr_recursive = "//*[local-name() = $name]"
name_expr_parent = ".//*[local-name() = $pname]/*[local-name() = $name]"

class resultHandler:
    """ Element result or resultSet element handling
    """
    def __init__(self, result_element):
        self._element = result_element
        self._result = ET.QName(result_element).localname == "result"
        
    def find(self, path):
        additional_elements = [c for ele in self._element.find(path) for c in ele.getchildren() if not ele.text.strip()]
        return pd.Series({ET.QName(k).localname: k.text 
                          for k in [*self._element.find(path), *additional_elements]
                          if k.text and k.text.strip()}).apply(pd.to_numeric, errors="ignore")
    @property
    def method(self):
        return self.find("./{*}propertyGroups/{*}common")
    
    @property
    def common_info(self):
        return self.find("./{*}common")
    
    @property
    def analysis_results(self):
        if not self._result:
            return pd.Series(dtype=object)
        d = {}
        for item in self._element.findall(".//{*}ApplicationItemReferenceInfo"):
            k, v = [c.text for c in item.findall(".//{*}NameValueMetaData/{*}Value")]
            d[k] = v
        return pd.Series(d).astype(int)
    
    @property
    def sample_info(self):
        if not self._result:
            return pd.Series(dtype=object)
        return self.find("./{*}sample/{*}propertyGroups/{*}common")
    
    @property
    def acquisition(self):
        if not self._result:
            return pd.Series(dtype=object)
        return self.find("./{*}sample/{*}propertyGroups/{*}acquisition")
    
    @property
    def references(self):
        return pd.DataFrame([dict(reference.attrib) for reference in self._element.find("./{*}references")])
   
    @property
    def results(self):
        return pd.concat([self.method, self.common_info, self.analysis_results, self.sample_info, self.acquisition])

class uepHandler:
    def __init__(self, file, index=None, mute=None):
        self._file = Path(file)
        self._index = index
        
        self._system_info = None
        self._sample_list = None
        self._streams = None
        self._streamLocation = {}
        
        self._peak_table = None
        self._content = None
    
        self._streams_relation = []
        self._section_relation = []
    
        self.mute = mute
    def __repr__(self):
        return f"{self.__class__.__name__}(r'{self._file}')"
        
    ### Properties
    @property
    def file(self):
        return self._file
    
    @property
    def index(self):
        if self._index is None:
            self._index = self._get_index()
        return self._index
    
    @property
    def content(self):
        if self._content is None:
            self._read_package_content()
        return self._content
    
    @property
    def streams_relation(self):
        self.content
        return 

    @property
    def section_relation(self):
        self.content
        return pd.DataFrame(self._section_relation)

    @property
    def streams(self):
        if self._streams is None:
            self.content
        return self._streams
            
    @property
    def file_count(self):
        return self.index["name"].str.split("_").apply(
                lambda x: ".".join(x[-1:])).value_counts().rename("file_count")
    
    @property
    def sample_list(self):
        if self._sample_list is None:
            self.content
            
        columns = ['file','name','systemName','sampleNo', 'sampleName', 'sampleId', 'acquisitionRunTime', 'acquisitionStartTime',
                   'submittedBy', 'replicateNumber', 'wellPosition',
                   'ecordId', 'eCordType', 'sampleType', 'sampleLevel',
                   'sampleWeight', 'dilution', 'injectionVolume', 'processingOptions',
                   'parentItemID', 'SubItemVersionNo', "nodeID", "itemID"
                   ]
        
        return self._sample_list.reindex(columns=columns)

    @property
    def system_info(self):
        if self._system_info is None:
            self.content
        temp = self.content.loc[lambda x: x["systemName"].notna(), self._system_info.index.difference(["category", "origin"])].drop_duplicates()
        assert len(temp) < 2, "Flere systemer?!"
        return self._system_info

    @property
    def references(self):
        lst = []
        for ix, row in self.content.dropna(subset=["references"]).iterrows():
            lst.append(row["references"].assign(index=ix))
        return pd.concat(lst)
    
    @property
    def parent_resultSet(self):
        resultSet_name = self.references.loc[lambda x: x["name"] == 'parentResultSetRef', "itemName"]
        if not resultSet_name.empty:
            return resultSet_name.iloc[0] if len(resultSet_name) == 1 else resultSet_name
    
    @property
    def peak_table(self):
        if self._peak_table is None:
            ms_cols = ["accurate",
                       "saturated",
                       "lockMassCorrected"]
            
            lst = []
            for (sampleid), streams in self.streams.loc[
                    lambda x: x["Name"].isin(['3D mass peak list', 'AMRTs from : 3D mass peak list'])
                    ].groupby("sampleId"):
                low_ms, high_ms, low_amrt, high_amrt = streams.index[-4:]
                
                for i, (ms_idx, amrt_idx) in enumerate(((low_ms, low_amrt),
                                                        (high_ms, high_amrt),
                                                        ), 1):
                    
                    MSPeakListStream = self.index_stream(ms_idx)
                    AMRTStream, peakStream, clusterStream = self.index_stream(amrt_idx)
                    
                    table = peakStream.join(MSPeakListStream.set_index('peakId')[ms_cols], 
                                             on='massPeakID')
                    table.insert(0, 'sampleId', sampleid)
                    table.insert(1, 'channel', i)
                    lst.append(table)
                
            if not lst:
                print(f"No peaks in {self._file}")
                return
            
            peak_table = pd.concat(lst, sort=False, ignore_index=True)
            
            assert all(peak_table.fragmentationMode == False), "fragmentationMode was not False"
            assert all(peak_table.isotopeID == 0), "isotopeID was not 0"
            assert all(peak_table.lockMassCorrected == True), "Not lockMassCorrected"
            
            self._peak_table = peak_table.drop([
                                    "isotopeID", 
                                    "driftTime", "driftTimeSD", "driftTimeInMilliSeconds", 
                                    "collisionCrossSection", "reducedCollisionCrossSection",
                                    "clusterLiftOffRT", "clusterTouchDownRT",
                                    "fragmentationMode",
                                    "lockMassCorrected",
                                    ], axis=1)      
        return self._peak_table
    

    def read_index(self, idx, parse_xml=True):
        with self.file.open("rb") as buffer:
            start, size = self.index.loc[idx, ["start", "size"]].values
            buffer.seek(start)
            bytestr = buffer.read(size)
            if parse_xml and self.index.loc[idx, "name"].endswith(".xml"):
                return ET.fromstring(bytestr)
            return bytestr
    
    def read_index_named(self, name, parse_xml=True, return_idx=False,):
        idx, = self.index.loc[lambda x: x["name"] == name].index
        
        if return_idx:
            return idx, self.read_index(idx)
        else:
            return self.read_index(idx)
    
    def read_section(self, itemSectionID):
        idx, = self.index.loc[lambda x: x["name"] == fr"\Data\{convert_str_id(itemSectionID)}_SectionXml.xml"].index
        return self.read_index(idx)
    
    def read_location(self, location):
        indices = self.index.loc[lambda x: (x["name"].str.endswith(location.upper() + "."))].index
        *other, idx = indices        
        
        if other and not self.mute:
            print(f"WARNING: {len(other)+1} matching locations!")
            
        return self.read_index(idx)
    
    ### Specialized methods
    def index_stream(self, idx):
        name, description, locations = self.streams.loc[idx, ["Name", "Description", "Location"]]
        lst = []

        for loc in locations:
            bytestream = self.read_location(loc)
            headers, fmt = self._fmt_dict[loc]
            try:
                df = pd.DataFrame(list(iter_unpack(fmt, bytestream)), columns=headers)
                lst.append(df)
            except Exception as e:
                print(f"Unpacking of {': '.join(self.streams.loc[idx, ['Sample', 'Name']].values)} failed!")
                raise e
            
            if description == "MS_SCANNING_BLOCKED":
                stream = self._streamLocation[loc]
                new_locs = [s.get("location").upper() for s in stream.xpath("./*[name() = $name]", name="stream")]
                for new_loc in new_locs:
                    if loc != new_loc:
                        return self._unpack_blocked(df, new_loc, loc)
        
        if name == "MSeBinning":
            t1, t2 = lst
            
        return lst if len(lst) > 1 else lst[0]
    
    ### Private methods   
    def _get_index(self):
        with self.file.open("rb") as buffer:
            buffer.seek(-250, 2)
            index_info = buffer.read()
            assert sum(list(index_info)) > 0, f'No UEP index: "{self.file}"'
            
            temp1, start, temp, size, *rest = index_info.split(b'"') #LastIndex, fileIndex comes right after.

            start = int(start) + int(size)
            buffer.seek(start)
            index = ET.fromstring(buffer.read()[:-174])
            return pd.DataFrame([dict(c.items()) for c in index]).astype({"start": np.int64, "size": np.int64})            
        
    def _read_package_index(self):
        package_index = self.read_index_named("\PackageIndex.xml")
        lst = []
        for entry in package_index.xpath(name_expr, name="indexEntry"):
            d = {**entry.attrib,}
            
            for child in entry:
                name = local_name(child)
                d[name] = child.text
            lst.append(d)
        
        return pd.DataFrame(lst)
        
    def _read_package_items(self, items):
        lst = []
        for item in items:
            lst.append(dict(item.attrib))
            for subversion in item.xpath(name_expr, name="subVersion"):
                lst.append(dict(subversion.attrib))
        data = (pd.DataFrame.from_dict(lst).rename(columns={"itemId": "itemID", "versionno": "versionNo"})
                                           .drop_duplicates("itemID", keep="last")
                                           .merge(self._read_package_index(), how="left")
                  )
        return data   
    
    def _read_itemMetaData(self, item_metadata):
        dct = {local_name(c): c for c in item_metadata}
        adminData = {local_name(c): c.text for c in dct.get("adminData", {}) if c.text.strip()}
        ID = {k: v for k, v in adminData.items() if "ID" in k}
        streams = [{**ID, **c.attrib} for c in dct.get("referencedStreams", {})]
        sections = [{**ID, **c.attrib} for c in dct.get("referencedSections", {})]
        if streams:  self._streams_relation.extend(streams)
        if sections: self._section_relation.extend(sections)
        return adminData
    
    def _read_package_entries(self, entries):
        entry_list = []
        subitem_list = []
        for entry in entries:
            instanceData, itemMetaData, subItems = entry
            entry_dict = {**{k.replace("Id", "ID"): v for k,v in entry.attrib.items()},
                          "instanceData": instanceData.text,
                          "itemMetaData": itemMetaData.text,
                          "subItems": len(subItems),
                          }
            metadata = self._read_itemMetaData(self.read_index_named(itemMetaData.text))
            entry_dict.update(metadata)
            if entry_dict["state"] == "SUPERCEDED":
                continue
            
            instance_data = self.read_index_named(instanceData.text)
            entry_dict["itemType"] = local_name(instance_data.tag)
            sample = instance_data.xpath(name_expr_parent, pname="sample", name="propertyGroups")
            if sample and local_name(instance_data.tag) == "result":
                sample, = sample
                series = pd.concat([pd.Series({local_name(c.tag): c.text for c in c}, dtype=object) for c in sample])
                entry_dict.update(series)
            
            elif sample and local_name(instance_data.tag) == "acquisitionList":
                pass
            
            elif local_name(instance_data.tag) == "resultSet":
                res = resultHandler(instance_data)
                self._system_info = res.results
                entry_dict.update(self._system_info)                 
                entry_dict["references"] = res.references
            
            subitems = pd.DataFrame([(entry_dict["itemID"], k.attrib.get("itemID"), k.attrib.get("versiono")) for k in subItems], 
                                     columns=["parentItemID", "itemID", "SubItemVersionNo"],
                                     )
            
            subitem_list.append(subitems)
            entry_list.append(entry_dict)

        subitems = pd.concat(subitem_list, sort=False, ignore_index=True).rename(columns={"nodeId":"nodeID"})
        data = (pd.DataFrame(entry_list).rename(columns={"itemId":"itemID"})
                                        .merge(subitems, how="left")
                                        )
        return data
    
    def _read_package_streams(self, streams):
        lst = []
        for stream in streams:
            d = {**stream.attrib}
            for child in stream:
                d[local_name(child)] = child.text
            meta = self.read_index_named(d.get('streamMetaDataExportPath'))
            d = {**d, **meta.attrib}
            lst.append(d)
    
        return pd.DataFrame(lst)
    
    def _read_section_metadata(self):
        """ Lettere at parse dem her, end InstanceData...
        """
        lst = []
        for i, row in self.index.loc[lambda x: x["name"].str.endswith("_SectionMetaData.xml")].iterrows():
            lst.append({c.tag: c.text for c in self.read_index(i)})
        return pd.DataFrame(lst)
    
    def _read_package_content(self):
        pc_element = self.read_index_named("\PackageContent.xml")
        comp, = pc_element.xpath(name_expr, name="compression")
        if not comp.text == "NONE":
            raise NotImplementedError("Reading compressed uep files is not implemented")
        
#        _package_items   = self._read_package_items(  pc_element.xpath(name_expr, name="item"))
        self._content = self._read_package_entries(pc_element.xpath(name_expr, name="entry"))
#        _package_streams = self._read_package_streams(pc_element.xpath(name_expr, name="stream"))
        
#        streams = _package_streams.merge(pd.DataFrame(self._streams_relation))
        sections = pd.DataFrame(self._section_relation).drop_duplicates().merge(self._read_section_metadata().rename(columns={"InternalSectionID": "internalSectionID"}))
        
        lst = []
        
        for i, row in sections.loc[lambda x: x["SectionType"].str.upper() == "CHANNELS"].iterrows():
            channels = self.read_section(row["internalSectionID"])
            info = self._build_stream_dict(channels).assign(**row)
            lst.append(info)
        
        self._sample_list = (self._content.loc[lambda x: (x["itemTypeID"] == "27fc931a-0b58-4632-855b-6ed7a1825d08")].copy()
                             .assign(**self.system_info[["name","systemName"]], file=str(self._file)))
            
        self._sample_list["sampleNo"] = range(1, len(self._sample_list)+1)
        self._streams = pd.concat(lst, ignore_index=True).merge(self._sample_list)

        self._fmt_dict = {k: get_headers_fmt(v) for k,v in self._streamLocation.items()}
        
        assert self._content["nodeID"].nunique() < 2, "Flere nodes?.."

    #TODO: Currently blind towards "<streams>"
    def _build_stream_dict(self, channels):
        data = []
        for channel in channels:
            name, = channel.xpath(name_expr, name="channelName")
            description = channel.get("dataDescription")
            
            locs = []
            if description == "MS_SCANNING_BLOCKED":
                streams = channel.xpath("./*[name() = $name]", name="stream")
            else:
                streams = channel.xpath(".//*[name() = $name]", name="stream")
                
            for stream in streams:
                loc = stream.get("location", "").upper()

                if not loc in self._streamLocation:
                    self._streamLocation[loc] = stream
                locs.append(loc)
                
            data.append((name.text, description, pd.Series(locs).unique().tolist()))
        return pd.DataFrame(data, columns=["Name", "Description", "Location"])
