import xml.etree.ElementTree as ET
import arcpy
class xml2csv:
	def __init__(self,xmlPath):
		self.fileAuthor = self.fileDOI = self.fileVersion_mj = self.fileVersion_mi  = self.ReSpecThVersion = self.firstPubDate = self.lastModDate = self.biblio_description = self.biblio_refDOI = self.biblio_loc = self.biblio_fig = self.bibliography_author = self.bibliography_journal = self.bibliography_pages = self.bibliography_title = self.bibliography_vol = self.bibliography_year = self.exp_type = self.appratus_kind = self.appratus_mode = self.pressure = self.P_unit = self.P_id = self.temperature = self.T_unit = self.T_id = self.response_type = self.response = self.response_unit = self.ignition_type = self.comment = self.data_point = None	
		self.xmlPath = xmlPath
		print(xmlPath)
		self.tree = ET.parse(self.xmlPath)
		self.root = self.tree.getroot()
		Target_kind = root.find(".//kind").text
		Target_mode = root.find(".//mode").text
		print(Target_kind)
		print(Target_mode)
		#for child in self.root:
			#if child.tag == "fileAuthor":
			#	self.fileAuthor = child.text
			#if child.tag == "fileDOI":
			#	self.fileDOI = child.text
			#if child.tag == "fileVersion":
				
			#if child.tag == "ReSpecThVersion":
			#if child.tag == "firstPublicationDate":
			#if child.tag == "lastModificationDate":
			#if child.tag == "bibliographyLink":
			#if child.tag == "experimentType":
			#if child.tag == "apparatus":
			#if child.tag == "commonProperties":
			#if child.tag == "dataGroup":
			#if child.tag == "ignitionType":
			#if child.tag == "comment":
			
			
            	
