from __future__ import print_function
import numpy
import matplotlib.pyplot as plt
import os
import string
from matplotlib.patches import Rectangle
import matplotlib.patches
import collections as c
import sys
import matplotlib.backend_bases as event_plot
from mpldatacursor import datacursor
from matplotlib.widgets import Slider as matplotSlider

#os.chdir("/Users/saljh8/Desktop/Code/AltAnalyze/Config/ExonViewFiles")

class EnsemblRegionClass:
    def __init__(self,start,stop,ensembl_exon,exon_region,strand):
        self.start = start
        self.stop = stop
        self.ensembl_exon = ensembl_exon
        self.exon_region = exon_region
        self.strand = strand
    def Start(self): return int(self.start)
    def Stop(self): return int(self.stop)
    def EnsemblExon(self): return self.ensembl_exon
    def EnsemblRegion(self): return self.exon_region
    def ExonBlock(self):
        return string.split(self.EnsemblRegion(),'.')[0]
    def Strand(self): return self.strand
    def Length(self): return abs(int(self.Start())-int(self.Stop()))
    def setChr(self, chr): self.chr = chr
    def Chr(self): return self.chr

class SplicingIndexClass:
    def __init__(self, reg_call, splicing_index, p_val, midas):
        self.reg_call = reg_call
        self.splicing_index = splicing_index
        self.p_val = p_val
        self.midas = midas
    def RegCall(self): return self.reg_call
    def SplicingIndex(self): return self.splicing_index
    def PVal(self): return self.p_val
    def Midas(self): return self.midas

def ProteinCentricIsoformView(Selected_Gene):
    Transcript_List = []
    Transcript_db = {}
    Exon_db = {}
    for line in open(Transcript_Annotations_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        if(line[0] == Selected_Gene):
            transcriptID = line[-1]
            exonID = line[5]
            start = line[3]
            stop = line[4]
            strand = line[2]
            chr = line[1]
            if 'chr' not in chr:
                chr = 'chr'+chr
            exon_data = EnsemblRegionClass(start,stop,exonID,None,strand)
            exon_data.setChr(chr)
            Transcript_List.append((transcriptID, exonID))
            try:
                Transcript_db[transcriptID].append(exon_data)
            except Exception:
                Transcript_db[transcriptID]=[exon_data]
            try:
                Exon_db[exonID].append(transcriptID)
            except Exception:
                Exon_db[exonID]=[transcriptID]
    
    Transcript_Protein_db = {}
    Protein_Transcript_db = {}
    Protein_List = []
    count = 0
    for line in open(Prt_Trans_File, "rU").xreadlines():
        if(count == 0):
            count = 1
            continue
        line = line.rstrip()
        line = line.split("\t")
        if(len(line) != 3):
            continue
        geneID = line[0]
        transcriptID = line[1]
        proteinID = line[2]
        if Selected_Gene == geneID:
            Transcript_Protein_db[transcriptID] = proteinID
            Protein_Transcript_db[proteinID] = transcriptID
            Protein_List.append(proteinID)
    
    Transcript_ExonRegion_db={}
    geneExonRegion_db={}
    exon_coord_db={}
    exonRegion_db={}
    AllBlocks = [("E", []), ("I", [])]
    # Store the exon region positions and later link them to the Ensembl exons
    for line in open(ExonRegion_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        geneID = line[0]
        exon_region = line[1]
        chr = line[2]
        strand = line[3]
        start = line[4]
        stop = line[5]
        er = EnsemblRegionClass(start,stop,exonID,exon_region,strand)
        if(geneID == Selected_Gene):
                Block_Num = exon_region[1:]
                I_E_id = exon_region[0]
                if(I_E_id == "E"):
                    AllBlocks[0][1].append(Block_Num)
                if(I_E_id == "I"):
                    AllBlocks[1][1].append(Block_Num)
                exon_added = False
                #Exon_List = line[7].split("|")
                exon_coord_db[chr,int(start),'start'] = exon_region
                exon_coord_db[chr,int(stop),'stop'] = exon_region
                exonRegion_db[Selected_Gene,exon_region] = er
                #print chr,start,'start'

    ETC_List = []
    for line in open(Etc_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        GeneLine = line[0].split(":")
        Gene = GeneLine[0]
        regcall = line[2]
        spl_index = line[3]
        pval = line[4]
        midas = line[5]
        S_I_data = SplicingIndexClass(regcall, spl_index, pval, midas)
        if(Gene == Selected_Gene):
            try:
                SplitGene = GeneLine[1].split("_")
                SplitGene = SplitGene[0]         
                ETC_List.append((SplitGene, S_I_data))
            except:
                pass

    ETC_dict = {}
       
    # Link the exon regions to the Ensembl exons
    for transcriptID in Transcript_db:
        for exon_data in Transcript_db[transcriptID]:
            start = exon_data.Start()
            stop = exon_data.Stop()
            chr = exon_data.Chr()
            strand = exon_data.Strand()
            start_exon_region = exon_coord_db[chr,start,'start']
            stop_exon_region = exon_coord_db[chr,stop,'stop']
            if '-' in strand:
                stop_exon_region,start_exon_region = start_exon_region,stop_exon_region
            regions = [start_exon_region]
            block,start_region = start_exon_region.split('.')
            start_region = int(float(start_region))
            block,stop_region = stop_exon_region.split('.')
            stop_region = int(float(stop_region))
            region = start_region+1
            while region<stop_region:
                er = block+'.'+str(region)
                regions.append(er)
                region+=1
            if stop_region != start_region:
                regions.append(stop_exon_region)
            for region in regions:
                er = exonRegion_db[Selected_Gene,region]
                try:
                    Transcript_ExonRegion_db[transcriptID].append(er)
                except:
                    Transcript_ExonRegion_db[transcriptID] = [er]
    
    exon_virtualToRealPos= c.OrderedDict()
    junction_transcript_db = {}
    for transcriptID in Transcript_ExonRegion_db:
            position=0
            Buffer=15
            for exon_object in Transcript_ExonRegion_db[transcriptID]:
                if position!=0:
                    if last_exon != exon_object.ExonBlock():
                        #print last_exon_region+'-'+exon_object.EnsemblRegion(),position,Buffer
                        junctionID = last_exon_region+'-'+exon_object.EnsemblRegion()
                        try: junction_transcript_db[transcriptID].append((position,position+Buffer, junctionID)) ### virtual junction positions
                        except: junction_transcript_db[transcriptID] = [(position,position+Buffer, junctionID)]
                        position+=Buffer
                        
                virtualStart = position
                virtualStop = virtualStart + exon_object.Length()
                position = virtualStop
                try:                    
                    exon_virtualToRealPos[transcriptID].append(([virtualStart,virtualStop],[exon_object.Start(), exon_object.Stop()],exon_object))
                except Exception:                    
                    exon_virtualToRealPos[transcriptID]=[([virtualStart,virtualStop],[exon_object.Start(), exon_object.Stop()],exon_object)]
                #print transcriptID,exon_object.ExonBlock(),exon_object.EnsemblExon(),exon_object.EnsemblRegion(),exon_object.Start(),exon_object.Stop(),virtualStart,virtualStop,"\n"
                last_exon = exon_object.ExonBlock()
                last_exon_region = exon_object.EnsemblRegion()

    for i in ETC_List:
        Region = i[0]
        S_I = i[1]
        Region = Region.split("-")
        if(len(Region) > 1):

            #Delete underscores from values.

            R_Start = Region[0]
            R_End = Region[1]
            R_Start = R_Start.split("_")
            R_End = R_End.split("_")
            R_Start = R_Start[0]
            R_End = R_End[0]
            R_Final = R_Start + "-" + R_End
            R_Type = R_Final[0]
            #print(R_Final)
            ETC_dict[R_Final] = S_I
            
        else:
            Region = Region[0]
            Region = Region.split("_")
            Region = Region[0]
            Region_type = Region[0]
            ETC_dict[Region] = S_I
            #if(Region_type == "E"):
            #    for entry in AllBlocks[0][1]:
            #        if(Region[1:] == entry):
            #            ETC_dict[("E" + entry)] = S_I
            #if(Region_type == "I"):
            #    for entry in AllBlocks[1][1]:    
            #        if(Region[1:] == entry):
            #            ETC_dict[("I" + entry)] = S_I

    #for a in ETC_dict:
    #        print(ETC_dict[a].RegCall(), a)
          
    #for i in junction_transcript_db:
    #    print i, junction_transcript_db[i], "\n"
    
    Protein_Pos_Db = {}
    last_protein=None
    stored_stop=None
    for line in open(Prt_Boundaries_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        proteinID = line[0]
        if(proteinID in Protein_List):
            Stop = int(line[-1])
            Start = int(line[-2])
            if(proteinID != last_protein):
                if stored_stop !=None:
                    #print proteinID,stored_start,stored_stop
                    Protein_Pos_Db[last_protein] = [[stored_start,stored_stop,None]]
                stored_start = int(Start)
            if(proteinID == last_protein):
                stored_stop = int(Stop)
            last_protein = str(proteinID)
    
    Protein_Pos_Db[last_protein] = [(stored_start,stored_stop,None)]
    Protein_virtualPos = RealToVirtual(Protein_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db)    
    
    domainAnnotation_db={}
    Domain_Pos_Db={}
    for line in open(Prt_Regions_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        proteinID = line[0]
        if proteinID in Protein_Pos_Db:
            domain_start = int(float(line[3]))
            domain_stop = int(float(line[4]))
            domainID = line[-2]
            domainName = line[-1]
            try:
                Domain_Pos_Db[proteinID].append((domain_start,domain_stop,domainID))
            except:
                Domain_Pos_Db[proteinID] = [(domain_start,domain_stop,domainID)]
            domainAnnotation_db[domainID] = domainName
    # Do the same for domain coordinates
    Domain_virtualPos = RealToVirtual(Domain_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db)
    return_val = ((junction_transcript_db, Protein_virtualPos, Domain_virtualPos, Transcript_db, exon_virtualToRealPos, ETC_dict))
    return return_val
    
def RealToVirtual(Protein_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db):
    Transcript_to_Protein_Coords = {}
    for proteinID in Protein_Pos_Db:
        transcript = Protein_Transcript_db[proteinID]
        strand = Transcript_ExonRegion_db[transcript][0].Strand()
        e_coords = exon_virtualToRealPos[transcript]
        if proteinID in Protein_Pos_Db:
            for q in Protein_Pos_Db[proteinID]:
                #print("Protein: ", proteinID)
                p_start = q[0]
                p_stop = q[1]
                annotation = q[2]
                if '-' not in strand:
                    p_start +=1
                    p_stop -= 1 ### seems to be off by 1
                virtual_p_start = None
                virtual_p_stop = None
                #print("E", len(e_coords))
                #print("Protein: ", proteinID)
                
                for i in range(len(e_coords)):
                    e_virtual_start, e_virtual_stop = e_coords[i][0]
                    #print("Sub-E", e_virtual_start, e_virtual_stop)
                    e_real_start,e_real_stop = e_coords[i][1]
                    #print e_real_start,e_real_stop
                    e = [e_real_start,e_real_stop]
                    e.sort()
                    p = [p_start,p_stop]
                    p.sort()
                    coord = e+p
                    coord.sort()
                    if (p_start<e[1] and p_start>e[0]) or p_start==e[1] or p_start==e[0]:
                        if '-' in strand:
                            offset = e_real_stop-p_start
                        else:
                            offset = p_start-e_real_start
                        virtual_p_start = offset+e_virtual_start
                        #print("Final_Val", proteinID, virtual_p_start)
                    if (p_stop<e[1] and p_stop>e[0]) or p_stop==e[1] or p_stop==e[0]:
                        if '-' in strand:
                            offset = e_real_stop-p_stop
                        else:
                            offset = p_stop-e_real_start
                        virtual_p_stop = offset+e_virtual_start
                if annotation != None:
                    #print("Entered", proteinID, virtual_p_start)
                    try:
                        Transcript_to_Protein_Coords[transcript].append((proteinID, annotation, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1]))
                    except Exception:
                        Transcript_to_Protein_Coords[transcript] = [(proteinID, annotation, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1])]
                else:
                    #print("Entered_2", proteinID, virtual_p_start)
                    Transcript_to_Protein_Coords[transcript] = proteinID, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1]
                #print transcript, proteinID, virtual_p_start, virtual_p_stop, p_start,p_stop, e_coords[0][0][0],e_coords[-1][0][1],annotation
    return Transcript_to_Protein_Coords
    
def remoteGene(gene):
    global Transcript_Annotations_File
    global ExonRegion_File
    global Selected_Gene
    global Prt_Trans_File
    global Prt_Regions_File
    global Prt_Boundaries_File
    global Etc_File
    import unique
    Selected_Gene = gene
    ExonRegion_File = unique.filepath("ExonViewFiles/Hs_Ensembl_exon.txt")
    Transcript_Annotations_File = unique.filepath("ExonViewFiles/Hs_Ensembl_transcript-annotations.txt")
    Prt_Trans_File = unique.filepath("ExonViewFiles/Hs_Ensembl_Protein__65_37.txt")
    Prt_Regions_File = unique.filepath("ExonViewFiles/Hs_ProteinFeatures_build_65_37.txt")
    Prt_Boundaries_File = unique.filepath("ExonViewFiles/Hs_ProteinCoordinates_build_65_37.tab")
    Etc_File = unique.filepath("ExonViewFiles/Hs_RNASeq_K562_SRSF2_P95mut_vs_K562_SRSF2_WT.ExpCutoff-5.0_average-splicing-index-ProcessedSpliceData.txt")
    #"ENSG00000005801"
    #"ENSG00000110514"
    total_val = ProteinCentricIsoformView(Selected_Gene)
    junctions = total_val[0]
    p_boundaries = total_val[1]
    p_domains = total_val[2]
    transcript_db = total_val[3]
    exon_db = total_val[4]
    splice_db = total_val[5]

    #for i in exon_db["ENST00000349238"]:
    #    print(i[2].EnsemblRegion())
    
    domain_color_list = []
    for i in p_domains:
        ploy = p_domains[i]
        for a in ploy:
            domain_color_list.append(a[1])

    domain_color_list = list(set(domain_color_list))
    domain_color_key = {}
    c_color1 = [0.8, 0.6, 0.1]
    c_color2 = [0.1, 0.6, 0.8]
    c_color3 = [0.6, 0.1, 0.8]
    c_color4 = [0.95, 0.6, 0.3]
    c_color5 = [0.3, 0.6, 0.95]
    c_color6 = [0.6, 0.3, 0.95]
    FLAG = 1

    for item in domain_color_list:
        if(FLAG == 1):
            domain_color_key[item] = c_color1
            FLAG = FLAG + 1
            continue
        if(FLAG == 2):
            domain_color_key[item] = c_color2
            FLAG = FLAG + 1
            continue
        if(FLAG == 3):
            domain_color_key[item] = c_color3
            FLAG = FLAG + 1
            continue
        if(FLAG == 4):
            domain_color_key[item] = c_color4
            FLAG = FLAG + 1
            continue
        if(FLAG == 5):
            domain_color_key[item] = c_color5
            FLAG = FLAG + 1
            continue
        if(FLAG == 6):
            domain_color_key[item] = c_color6
            FLAG = 1
            continue

    for i in domain_color_key:
        print(i, domain_color_key[i], "\n")
    
    Y = 50
    Transcript_to_Y = {}
    for transcript in transcript_db:
        Transcript_to_Y[transcript] = Y
        Y = Y + 200
    import traceback

    ylim = Y + 200
    currentAxis = plt.gca()
    ax = plt.axes()
    X_Pos_List = []
    for transcript in transcript_db:
        try:
            Junc_List = junctions[transcript]
            y_pos = Transcript_to_Y[transcript]
            Gene_List = exon_db[transcript]
            color_flag = 1
            for entry in Gene_List:
                G_start = entry[0][0]
                G_end = entry[0][1]
                Exon_Object = entry[2]
                try:
                    LabelClass = splice_db[Exon_Object.EnsemblRegion()]
                    ExonName = Exon_Object.EnsemblExon()
                    RegCall = LabelClass.RegCall()
                    SplicingIndex = LabelClass.SplicingIndex()
                    PVal = LabelClass.PVal()
                    Midas = LabelClass.Midas()
                    Label = "\n" + "Exon: " + str(ExonName) + "\n" + "RegCall: "  + str(RegCall) + "\n" + "Splicing Index: " + str(SplicingIndex) + "\n" + "P-Value: " + str(PVal) + "\n" + "Midas Value: " + str(Midas)
                    if(RegCall == "UC"):
                        color_choice = "Grey"
                    else:
                        S_Int = float(SplicingIndex)
                        if(S_Int > 0):
                            color_choice = (0.7, 0.7, 0.99)
                        if(S_Int < 0):
                            color_choice = (0.8, 0.4, 0.4)
                                            
                except:
                    #print(traceback.format_exc());sys.exit()
                    Label = ""
                    color_choice = "Grey"
                #print("Start", G_start, "end", G_end, "Region", entry[2].EnsemblRegion())
                if((color_flag % 2) == 0):
                    currentAxis.add_patch(Rectangle((G_start, y_pos), (G_end - G_start), 50, color = color_choice, label = (entry[2].EnsemblRegion() + Label)))
                if((color_flag % 2) != 0):                   
                    currentAxis.add_patch(Rectangle((G_start, y_pos), (G_end - G_start), 50, color = color_choice, label = (entry[2].EnsemblRegion() + Label)))
                color_flag = color_flag + 1
                
                

            for entry in Junc_List:
                try:
                    LabelClass = splice_db[entry[2]]
                    RegCall = LabelClass.RegCall()
                    SplicingIndex = LabelClass.SplicingIndex()
                    PVal = LabelClass.PVal()
                    Midas = LabelClass.Midas()
                    Label = "\n" + "RegCall: " + str(RegCall) + "\n" + "Splicing Index: " + str(SplicingIndex) + "\n" + "P-Value: " + str(PVal) + "\n" + "Midas Value: " + str(Midas)
                    if(float(SplicingIndex) > 0):
                        color_junc = "blue"
                    if(float(SplicingIndex) < 0):
                        color_junc = "red"
                    if(RegCall == "UC"):
                        color_junc = "grey"
                except:
                    Label = ""
                    color_junc = "grey"
                currentAxis.add_patch(Rectangle((entry[0], y_pos), (entry[1] - entry[0]), 50, color = "White", label = (str(entry[2]) + Label)))
                ax.arrow(entry[0], (y_pos+50), 8, 40, label = (str(entry[2]) + Label), color = color_junc)
                ax.arrow((entry[0] + 8), (y_pos+90), 11, -40, label = (str(entry[2]) + Label), color = color_junc)


            P_Bound_List = p_boundaries[transcript]
            P_Domain_List = p_domains[transcript]
            E_Start = P_Bound_List[-2]
            E_End = P_Bound_List[-1]
            P_Start = P_Bound_List[1]
            P_End = P_Bound_List[2]
            #print("Boundaries: ", P_Start, P_End)
            X_Pos_List.append(int(E_End))
            #currentAxis.add_patch(Rectangle((E_Start, y_pos), E_End, 50, color = "Blue"))
            try:
                currentAxis.add_patch(Rectangle((P_Start, (y_pos + 120)), (P_End - P_Start), 10, label = ("Protein: " + str(P_Bound_List[0]))))
            except:
                pass
            for entry in P_Domain_List:
                #print("Domain", entry)
                color_domain_choice = domain_color_key[entry[1]]
                currentAxis.add_patch(Rectangle((entry[2], y_pos + 100), (entry[3] - entry[2]), 50, color = color_domain_choice, label= ("Protein: " + str(entry[0]) + "\n" + "Domain: " + str(entry[1]))))
        except:
            continue
    plt.ylim([0.0, ylim])
    try:
        max_x = max(X_Pos_List)
    except:
        max_x = 5000
    try:
        plt.xlim([0.0, max_x])
    except:
        plt.xlim([0.0, 3000])
    datacursor(hover=True, formatter='{label}'.format, bbox=dict(fc='yellow', alpha=1), arrowprops=None)
    plt.show()

        #currentAxis.add_patch(Rectangle

if __name__ == "__main__":
    #Selected_Gene = sys.argv[1]
    Selected_Gene = 'ENSG00000005801'
    remoteGene(Selected_Gene)