#!/usr/bin/python
from Bio.SeqUtils import GC
#import pylab
import os, time, subprocess
#import cProfile


def blat_against_human(TargetSeqs, Outfile, outType):
    os.environ["PATH"] = os.environ["PATH"] + ":/home/kaymazy/biotools/BLAT/"
    
    parameters = ['blat',
               '/home/kaymazy/genome.fa', #target Sequence
               '%s' %TargetSeqs, # Query Sequence
               '%s' %Outfile, # Output File
               '-maxIntron=2',
               '-out=%s' %outType] # output type
    #Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    return


def filter_blat_results(Blast8File, e_val_th, score_th, filtered_out_filename):
    
    filtered_out_file = open('%s' %filtered_out_filename, "w")
    
    with open("%s" %Blast8File, "r") as f_blast8:
        
        Blast8_lines = f_blast8.readlines()
        blast_hits_list = []
        for Blast8_line in Blast8_lines:
            q_name = str(Blast8_line.strip().split("\t")[0])
            hit_chr = str(Blast8_line.strip().split("\t")[1])
            q_str = str(Blast8_line.strip().split("\t")[6])
            q_end = str(Blast8_line.strip().split("\t")[7])
            hit_str = str(Blast8_line.strip().split("\t")[8])
            hit_end = str(Blast8_line.strip().split("\t")[9])
            e_val = float(Blast8_line.strip().split("\t")[10])
            score = float(Blast8_line.strip().split("\t")[11])
            
            if e_val < e_val_th and score > score_th:
                
                print e_val, score
                filtered_out_file.write(q_name+"\t"+q_str+"\t"+q_end+"\t"+hit_chr+"\t"+hit_str+"\t"+hit_end+"\t"+str(score)+"\t"+str(e_val)+"\n")
                blast_hits_list.append(q_name)
                
    filtered_out_file.close()       
    return blast_hits_list


def Linearize_genome(GenomeFasta_name):
    
    genomeLin = ''
    with open('%s' %GenomeFasta_name, "r") as GenomeFastaFile:
        FastaLines = GenomeFastaFile.readlines()
        for FastaLine in FastaLines:
            if FastaLine.startswith(">"):
                pass
            else:
                genomeLin = genomeLin + FastaLine.strip()
    return genomeLin


def get_HMRs(ProbeGenomeFasta_name, HitsFasta_filename, OutPsl, e_val_th, score_th, OutBlast8_filename, filtered_out_filename, hmr_bed_filename):
    
    Hitsoutfasta = open('%s' %HitsFasta_filename, "w")
    
    outType = 'psl'
    blat_against_human(ProbeGenomeFasta_name, OutPsl, outType)
    
    wt_genomeLin = Linearize_genome(ProbeGenomeFasta_name)
    
    hitDic = {}
    with open("%s" %OutPsl, "r") as f_outPsl:
        Psl_lines = f_outPsl.readlines()
        
        hitId = 0
        for Psl_line in Psl_lines:
            
            if Psl_line.startswith("psL") or\
             Psl_line.startswith("match") or\
              Psl_line.startswith("\n") or\
               Psl_line.startswith("\t") or\
                Psl_line.startswith("-") or\
                Psl_line.startswith(" "):
                pass
            else:
                match = Psl_line.strip().split("\t")[0]
                block_size = Psl_line.strip().split("\t")[18]
                blocks_count = Psl_line.strip().split("\t")[17]
                q_name = Psl_line.strip().split("\t")[9]
                q_start = int(Psl_line.strip().split("\t")[11])
                q_end = int(Psl_line.strip().split("\t")[12])
                hitId = hitId +1
                hitDic[hitId] = q_name, q_start, q_end
                for item in block_size.split(","):
                    if item == '':
                        pass
                    elif int(item) > 40:
                
                        Hitsoutfasta.write('>%d' %hitId+"\n"+wt_genomeLin[q_start:q_end]+"\n")
                        print hitId, item, block_size, blocks_count, match, q_start, q_end, wt_genomeLin[q_start:q_end]
                    else:
                        pass
                    
    Hitsoutfasta.close()
    
    outType = 'blast8'
    blat_against_human(HitsFasta_filename, OutBlast8_filename, outType)
    
    blast_hit_ids = filter_blat_results(OutBlast8_filename, e_val_th, score_th, filtered_out_filename)
    
    hmr_bed = open('%s' %hmr_bed_filename, "w")
    for item in hitDic:
        item = str(item)
        if item in blast_hit_ids:
            item = int(item)
            hmr_bed.write(str(hitDic[item][0])+"\t"+str(hitDic[item][1])+"\t"+str(hitDic[item][2])+"\n" )
            
    
    hmr_bed.close()
    f_outPsl.close()
    
    return


def design_probes(ProbeGenomeFasta_name, hmr_bed_filename, probes_txt_out_filename, probes_fa_out_filename, FirstGenomeID, X):

    probes_txt_out = open('%s' %probes_txt_out_filename, "w")
    probes_fa_out = open('%s' %probes_fa_out_filename, "w")
    
    genome = Linearize_genome(ProbeGenomeFasta_name)
    
    banned_coords_list = []
    with open("%s" %hmr_bed_filename, "r") as hmr_bed_file:
        linesH = hmr_bed_file.readlines()
        print " Yes, we can read hmr bed file"
        for lineH in linesH:
            hmr_start = int(lineH.strip().split("\t")[1])
            hmr_end = int(lineH.strip().split("\t")[2])
            for b in range(hmr_start-90, hmr_end-30, 1):
                banned_coords_list.append(b)
    
    banned_coords_set = set(banned_coords_list)
    
    probes_seq_list = []
    gc_list = []
    skipped_probe =0  
    p_id = 0
    m = 120/int(X)
    for n in range(0, 120, m):
        for i in range(n, len(genome), 120):
            #First tile set
            probe_seq = genome[i:i+120]
            count = 1
            
            if i in banned_coords_set:
                skipped_probe = skipped_probe +1
                #print probe_seq
                pass
            
            elif len(probe_seq) < 120:

                probe_seq = genome[len(genome)-120:len(genome)]
                
                probes_seq_list.append(probe_seq)
                
                p_id = int(p_id) + 1
                p_id = '%(p_id)012d' %{"p_id":p_id}
                probe_coor = '%(ID)s-%(str)d-%(end)d' %{"ID":FirstGenomeID,"str":len(genome)-120, "end":len(genome)}
                
            
                gc_content = GC(probe_seq) #Calculate GC content
                gc_list.append(gc_content)
            
                #Write probes into file format for SureDesign
                probes_txt_out.write(str(probe_coor)+"\t"+str(p_id)+"\t"+str(probe_seq)+"\t"+str(count)+"\t"+str(gc_content)+"\n")
                #probes_txt_out.write(str(probe_coor)+"\t"+str(p_id)+"\t"+str(probe_seq)+"\t"+str(count)+"\n")
                probes_fa_out.write(">"+str(probe_coor)+"\n"+str(probe_seq)+"\n")
            
                if probe_seq in probes_seq_list:
                    count = count + 1
                else:
                    pass
            else:
            
                probes_seq_list.append(probe_seq)
                probe_coor = '%(ID)s-%(str)d-%(end)d' %{"ID":FirstGenomeID,"str":i, "end":i+len(probe_seq)}

                p_id = int(p_id) + 1                
                p_id = '%(p_id)012d' %{"p_id":p_id}
            
                gc_content = GC(probe_seq) #Calculate GC content
                gc_list.append(gc_content)
                
                #Write probes into file format for SureDesign
                probes_txt_out.write(str(probe_coor)+"\t"+str(p_id)+"\t"+str(probe_seq)+"\t"+str(count)+"\t"+str(gc_content)+"\n")
                probes_fa_out.write(">"+str(probe_coor)+"\n"+str(probe_seq)+"\n")
            
                if probe_seq in probes_seq_list:
                    count = count + 1
                else:
                    pass
    
                
    probes_fa_out.close()
    probes_txt_out.close()
    print "Lenght of the Probe Genome", len(genome) 
    print "Number of Probes:", len(probes_seq_list)
    print "Number of Skipped Probes:", skipped_probe
    return gc_list 


def build_bwt_index(TargetGenomeFasta_name, Index_base):
    os.environ["PATH"] = os.environ["PATH"] + ":/home/kaymazy/pipeline/Bowtie/bowtie2-2.0.6/"
    parameters = ['bowtie2-build',
                  '%s' % TargetGenomeFasta_name,
                  '%s' % Index_base] #Index base name and directory
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    return


def bowtie_single(Index_base, TargetGenomeFasta_name, ProbesFasta, Out_sam_name):
    os.environ["PATH"] = os.environ["PATH"] + ":/home/kaymazy/pipeline/Bowtie/bowtie2-2.0.6/"
    os.environ["PATH"] = os.environ["PATH"] + ":/home/kaymazy/pipeline/SamTools/samtools-0.1.18/"
    
    parameters = ['bowtie2',
               '-f',
               '--end-to-end',
               '-x',
               '%s' %Index_base,
               '-U',
               '%s' %ProbesFasta, # Query Sequence
               '-S', # Output File
               '%s' %Out_sam_name]
    
    #Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'view', '-bS', '%s' %Out_sam_name, '>', '%s' %Out_sam_name[:-4]+'.bam']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'sort', '%s' %Out_sam_name[:-4]+'.bam', '%s' %Out_sam_name[:-4]+'_sorted']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'view', '-h', '%s' %Out_sam_name[:-4]+'_sorted.bam', '>', '%s' %Out_sam_name[:-4]+'_sorted.sam'] 
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'index', '%s' %Out_sam_name[:-4]+'_sorted.bam']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'fillmd', '-e','%s' %Out_sam_name[:-4]+'_sorted.bam', '%s' %TargetGenomeFasta_name, '>', '%s' %Out_sam_name[:-4]+'_MDfilled.sam']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'view', '-bS', '%s' %Out_sam_name[:-4]+'_MDfilled.sam', '>', '%s' %Out_sam_name[:-4]+'_MDfilled.bam']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'sort', '%s' %Out_sam_name[:-4]+'_MDfilled.bam', '%s' %Out_sam_name[:-4]+'_sorted_MDfilled'] 
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'index', '%s' %Out_sam_name[:-4]+'_sorted_MDfilled.bam']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    return


def parse_cigar_for_indel(cigar):
    """
    Parse CIGAR string and return pysam-suitable tuples.

    MIDNSHP => 0123456

    >>> parse_cigar("4S17M8D4M9I3H")
    [(4, 4), (0, 17), (2, 8), (0, 4), (1, 9), (5, 3)]
    """
    # maps operations to their integer encodings
    ops = dict( (i,j) for (j, i) in enumerate('MIDNSHP') )
    result = []
    n = ''
    
    for c in cigar:
        if c.isdigit():
            n += c
        elif c in ops:
            result.append( [ops[c], int(n)] )
            n = ''
    return result


def Probe_aln_filter(ProbesTxtFile, MDfilled_sam_file, InDelThreshold, PassedProbes_fa_out, Out_sam_name, SamForTilingDepth, TotalIndelMiss, MaxPerMatch):
    
    PassedSamFile = open('%s' %SamForTilingDepth, "w")
    PASSED_probes_IDs = []
    #Put original Probe sequences into a dictionary with their original unique IDs.
    sequences = {}
    
    with open("%s" %ProbesTxtFile, "r") as Txt:
        linesTxt = Txt.readlines()
        for lineTxt in linesTxt:
            ProbeId = str(lineTxt.strip().split("\t")[0])
            seq = lineTxt.strip().split("\t")[2]
            sequences[ProbeId] = ProbeId, seq

    ProbesPassed = open("%s" %PassedProbes_fa_out, "w")
    #Evaluate each probe according to total number of mismatch and indels
    Num_Probes_passed = 0
    #Open alignment file as SAM format
    with open("%s" %MDfilled_sam_file, "r") as f_MDFillsam:
        MDlines = f_MDFillsam.readlines()
        for MDline in MDlines:
            #Skip the header lines
            if MDline.startswith("@"):
                pass
            else:
                #store Reference sequence NAME
                Rname = str(MDline.strip().split("\t")[2])
                #Store "Modified" Probe sequence in a string variable. Modified sequence: Match bases are represented as "=".
                probeSeq = MDline.strip().split("\t")[9]
                #Store probe ID
                probeID = str(MDline.strip().split("\t")[0])
                #Count the total number of mismatches per probe sequence. "miss" is equal to total # of mismatches.
                miss = 120 - probeSeq.count("=")
                #Start counting indels per probe sequence.
                indel_count = 0
                InDel = [1, 2]
                #Retrieve the CIGAR string from the sam file.
                cigar = MDline.strip().split("\t")[5]
                #Parse CIGAR string and store in a list per probe sequence.                
                cigar_digit = parse_cigar_for_indel(cigar)
                #Count the total number of indels.
                #Take only aligned probes (without * tag) into count
                if Rname == "*":
                    pass
                else:
                    Bases = ['A','T','G','C']
                    start = 0
                    inserted_bases = 0
                    L_list = []
                    for i in range(len(cigar_digit)):
                        #Count only if Insertion or Deletion
                        if cigar_digit[i][0] in InDel:
                            #increase the indel # through out the probe sequence
                            indel_count = cigar_digit[i][1] + indel_count
                        L=0
                        if cigar_digit[i][0] == 1:
                            inserted_bases = inserted_bases + cigar_digit[i][1]
                        if cigar_digit[i][0] == 0:       
                            start = start + inserted_bases
                            end = start + cigar_digit[i][1]
                            for j in range(0, i+1):
                                for l in range(len(probeSeq[start:end])):
                                    if probeSeq[start:end][l] == "=":
                                        L = L+1
                                        if l == (len(probeSeq[start:end])-1):
                                            L_list.append(L)
                                            L=0
                                    elif probeSeq[start:end][l] in Bases:
                                        L_list.append(L)
                                        L = 0
                            start = start + cigar_digit[j][1] + inserted_bases
                            end = start + cigar_digit[i][1]
                    Lmax = max(L_list)
                    
                    #Filter probes using indel and mismatch thresholds
                    if (indel_count < InDelThreshold and miss < TotalIndelMiss) or MaxPerMatch < Lmax:
                    #if (indel_count < InDelThreshold and miss < MisMatchThreshold and (indel_count + miss) < TotalIndelMiss) or MaxPerMatch < Lmax:
                        #write probes' actual sequences as a fasta format output by retrieving from "sequences" dictionary.
                        ProbesPassed.write('>'+str(probeID)+"\n"+str(sequences[str(probeID)][1])+"\n")
                        #Store passed probe ids in a list
                        PASSED_probes_IDs.append(probeID)
                        #Count the total # of probes Passed the filtering.
                        Num_Probes_passed = Num_Probes_passed + 1
                        
                    #    if miss > TotalIndelMiss:
             #           print "Max L:", Lmax, probeSeq, "Total mismatch:", miss, "CIGAR string:", cigar_digit#, "Total InDelMis:", (indel_count + miss) 
                    else:
                        pass
        PASSED_probes_ID_set = set(PASSED_probes_IDs)
        
    
    with open("%s" %Out_sam_name[:-4]+'_sorted.sam', "r") as SamFile:
        Samlines = SamFile.readlines()
        for Samline in Samlines:
            samlineID = Samline.strip().split("\t")[0]
            
            if Samline.startswith("@"):
                PassedSamFile.write(Samline)
            elif samlineID in PASSED_probes_ID_set:
                PassedSamFile.write(Samline)
            else:
                pass   
        
    #PassedSamFile.close()
    #Out_sam_name.close()
    return Num_Probes_passed




def samtools_depth(TargetGenomeFasta_name, SamForTilingDepth, TilingDepth_out_filename, Gap_filled_coverage_out_filename):
    
    os.environ["PATH"] = os.environ["PATH"] + ":/home/kaymazy/pipeline/SamTools/samtools-0.1.18/"
    
    parameters = ['samtools', 'view', '-bS', '%s' %SamForTilingDepth,'>','%s' %SamForTilingDepth[:-4]+'.bam']               
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['samtools', 'depth', '%s' %SamForTilingDepth[:-4]+'.bam', '>', '%s' %TilingDepth_out_filename]               
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
        
    Targenome = Linearize_genome(TargetGenomeFasta_name)
    
    Gap_filled_coverage_out = open('%s' %Gap_filled_coverage_out_filename, "w")
    
    depth = {}
    coords = []
    with open("%s" %TilingDepth_out_filename, "r") as f_Depth:
        Depthlines = f_Depth.readlines()
        for Depthline in Depthlines:
            base_coord = Depthline.strip().split("\t")[1]
            coords.append(int(base_coord))
            tiling_depth = Depthline.strip().split("\t")[2]
            depth[int(base_coord)] = base_coord, tiling_depth
            
        coords_set = set(coords)
        for i in range(1, len(Targenome),1):   
            if i in coords_set:
                
                Gap_filled_coverage_out.write(str(depth[i][0])+"\t"+str(depth[i][1])+"\n")
                               
            else:
                
                Gap_filled_coverage_out.write(str(i)+"\t"+str(0)+"\n")
        
    Gap_filled_coverage_out.close()
    
    return


def new_probes_design(TargetGenomeFasta_name, hmr_bedForTarGet_filename, Gap_filled_coverage_out_filename, New_probes_txt_out_filename, New_probes_fa_out_filename, ID, X):
    
    New_probes_fa_out = open('%s' %New_probes_fa_out_filename, "w")
    New_probes_txt_out = open('%s' %New_probes_txt_out_filename, "w")
    
    base_list = []
    depth = []
    
    TargenomeLin = Linearize_genome(TargetGenomeFasta_name)
    print "Length of the target genome new_probe_design", len(TargenomeLin)
    
    Tarbanned_coords_list = []
    with open("%s" %hmr_bedForTarGet_filename, "r") as hmr_bedTar_file:
        linesTarH = hmr_bedTar_file.readlines()
        for lineTarH in linesTarH:
            tarhmr_start = int(lineTarH.strip().split("\t")[1])
            tarhmr_end = int(lineTarH.strip().split("\t")[2])
            for b in range(tarhmr_start-90, tarhmr_end-30, 1):
                Tarbanned_coords_list.append(b)
    
    Tarbanned_coords_set = set(Tarbanned_coords_list)
    
    
    with open("%s" %Gap_filled_coverage_out_filename,"r") as f_GapFilled:
        GapFlines = f_GapFilled.readlines()
        for GapFline in GapFlines:
            tiling_depth = GapFline.strip().split("\t")[1]
            base = GapFline.strip().split("\t")[0]
            base_list.append(base)
            depth.append(tiling_depth)
            
    new_probes_seq_list = []
    count_new_probes = 0 
    p_id = 0
    m = 120/int(X)
    for i in range(0, len(base_list), m):
        count = 0
        new_probe_seq = TargenomeLin[i:i+120]
        new_gc_content = GC(new_probe_seq)
        if i in Tarbanned_coords_set:
                
                print new_probe_seq
                pass
            
        elif depth[i] < str(2):
            new_probe_seq = TargenomeLin[i:i+120]
            new_gc_content = GC(new_probe_seq)
            if len(new_probe_seq) < 120:
                new_probe_seq = TargenomeLin[len(TargenomeLin)-120:len(TargenomeLin)]
                new_gc_content = GC(new_probe_seq)
                new_probes_seq_list.append(new_probe_seq)
                
                probe_coor = '%(ID)s-%(str)d-%(end)d' %{"ID":ID,"str":len(TargenomeLin)-120, "end":len(TargenomeLin)}
                
                p_id = int(p_id) + 1 + len(TargenomeLin)
                p_id = '%(p_id)012d' %{"p_id":p_id}
                
                count_new_probes = count_new_probes +1
                if new_probe_seq in new_probes_seq_list:
                    count = count + 1
                else:
                    pass
            else:
                probe_coor = '%(ID)s-%(str)d-%(end)d' %{"ID":ID,"str":i, "end":i+len(new_probe_seq)}
                new_probes_seq_list.append(new_probe_seq)                
                p_id = int(p_id) + 1 + len(TargenomeLin)
                p_id = '%(p_id)012d' %{"p_id":p_id}
                
                count_new_probes = count_new_probes +1
                if new_probe_seq in new_probes_seq_list:
                    count = count + 1
                else:
                    pass
            #Write probes into file format for SureDesign
            New_probes_txt_out.write(str(probe_coor)+"\t"+str(p_id)+"\t"+str(new_probe_seq)+"\t"+str(count)+"\t"+str(new_gc_content)+"\n")
            #Write probes into file  FASTA format
            New_probes_fa_out.write(">"+str(probe_coor)+"\n"+str(new_probe_seq)+"\n")
            
        else:
            pass
    New_probes_fa_out.close()
    New_probes_txt_out.close()
    return count_new_probes

def update_probeSet(ProbeSet, New_ProbeSet, Updated_ProbeSet):
    
    parameters = ['mv', '%s'%ProbeSet, '%s' %ProbeSet+'_2']
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    parameters = ['cat', '%s'%ProbeSet+'_2', '%s' %New_ProbeSet,'>','%s' %Updated_ProbeSet]               
    cmd = ' '.join(parameters)
    start = time.time()
    p = subprocess.Popen(cmd, shell=True)
    pid,ecode = os.waitpid(p.pid, 0)
    print round(time.time() - start,1)
    
    return

def Probes_count(Updated_ProbeSetTXT, probe_counts):
    
    Genomes = ['wt', 'AG876', 'Mutu', 'Akata', 'GD1']
    for item in Genomes:
        parameters = ['grep', '%s' %item, '%s' %Updated_ProbeSetTXT, '|', 'wc', '-l', '>>', '%s' %probe_counts ]               
        cmd = ' '.join(parameters)
        start = time.time()
        p = subprocess.Popen(cmd, shell=True)
        pid,ecode = os.waitpid(p.pid, 0)
        print round(time.time() - start,1)
    
    return


def boost(Updated_ProbeSetTXT, b_type, Pool_factor, Boosted_probes_txt_out_filename):
    
    BoostedProbes_txt_out = open('%s' %Boosted_probes_txt_out_filename, "w")
    #read the updated probe set
    factors = [2, 3, 8, 16, 4]
    with open('%s' %Updated_ProbeSetTXT, "r") as u_file:
        u_lines = u_file.readlines()
        for u_line in u_lines:
            gc = float(u_line.strip().split("\t")[4])
            seq = u_line.strip().split("\t")[2]
            name = u_line.strip().split("\t")[0]
            uid = u_line.strip().split("\t")[1]
            count = int(u_line.strip().split("\t")[3])
            
            count = Pool_factor * count
            
            if b_type == 'Maximum':
                if gc <= 59.00:#Sure
                    #print name, "\t", uid, "\t", seq, "\t", count, "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count)+"\n")
                    
                elif 62.00 >= gc > 59.00:#Sure
                    #print name, "\t", uid, "\t", seq, "\t", factors[0], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[0])+"\n")    
                
                elif 66.00 >= gc > 62.00:#Sure
                    #print name, "\t", uid, "\t", seq, "\t", factors[1], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[1])+"\n")
                
                elif 72.00 >= gc > 66.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[2], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[2])+"\n")
                
                elif gc > 72.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[3], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[3])+"\n")
                
                else:
                    pass
                
            elif b_type == 'Balanced':
                if gc < 65.00:
                    #print name, "\t", uid, "\t", seq, "\t", count, "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count)+"\n")
                    
                elif gc >= 65.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[0], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[0])+"\n")
                
            elif b_type == 'Mild_Balanced':
                if gc < 65.00:
                    #print name, "\t", uid, "\t", seq, "\t", count, "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count)+"\n")
                    
                elif 75.00 >= gc > 65.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[0], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[0])+"\n")
                
                elif 85.00 >= gc > 75.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[3], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[1])+"\n")
                    
                elif gc > 85.00:
                    #print name, "\t", uid, "\t", seq, "\t", factors[3], "\t", gc, "\n"
                    BoostedProbes_txt_out.write(name+"\t"+uid+"\t"+seq+"\t"+str(count * factors[4])+"\n")
                    
                else:
                    pass
            
    return
    

dirInput = '/home/kaymazy/EBV_Capture/EBV_Genomes/'
dirSubOuts = '/home/kaymazy/EBV_Capture/test1001/'
FirstGenomeID = 'wt'
TargetID = 'AG876'

InDelThreshold = 4
#MisMatchThreshold = 5 ### No need to specify ---> mismatches are included in indel count
TotalIndelMiss = 8
MaxPerMatch = 80
X = 4 #Tiling coverage
Boost_Type = 'Mild_Balanced'  #'Maximum' or 'Balanced'
Pool_factor = 6 

#FILES
# Design Probes based on Genomic Fasta file.
ProbeGenomeFasta_name = '%(dir)s%(FirstG)s_ebv.fa' %{"dir":dirInput, "FirstG":FirstGenomeID}
TargetGenomeFasta_name = '%(dir)s%(TarID)s_ebv.fa' %{"dir":dirInput, "TarID":TargetID}
OutPsl = '%(dir)sblat_%(FirstG)s_to_human.psl' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
HitsFasta_filename = '%(dir)sblat_%(FirstG)s_to_human_BlockSizeFiltered_HitSeqs.fasta' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
TargetSeqs = '%(dir)sblat_%(FirstG)s_to_human_BlockSizeFiltered_HitSeqs.fasta'%{"dir":dirSubOuts, "FirstG":FirstGenomeID}
OutBlast8_filename = '%(dir)sblat_%(FirstG)s_to_human_BlockSizeFiltered.blast8' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
filtered_out_filename = '%(dir)sblat_%(FirstG)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_blast8.bed'%{"dir":dirSubOuts, "FirstG":FirstGenomeID}
hmr_bed_filename = '%(dir)sblat_%(FirstG)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_HMRs.bed' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
#Produced Probes out file in SureDesign format
probes_txt_out_filename = '%(dir)s%(FirstG)s_Probes.txt' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
#Produced Probes out file in Fasta format
probes_fa_out_filename = '%(dir)s%(FirstG)s_Probes.fasta' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}

#Step 1.
#get HMR
################################cProfile.run('get_HMRs(ProbeGenomeFasta_name, HitsFasta_filename, OutPsl, 0.0001, 40.0, OutBlast8_filename, filtered_out_filename, hmr_bed_filename)')
#design Probes
design_probes(ProbeGenomeFasta_name, hmr_bed_filename, probes_txt_out_filename, probes_fa_out_filename, FirstGenomeID, X)# X=2 for 2x coverage
#return ProbeSet


Index_base = '%(dir)s%(TarID)s_ebv' %{"dir":dirSubOuts, "TarID":TargetID}
ProbesFasta = '%(dir)s%(FirstG)s_Probes.fasta' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
Out_sam_name = '%(dir)s%(FirstG)sProbes_aligned_to_%(TarID)s.sam' %{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
ProbesTxtFile = '%(dir)s%(FirstG)s_Probes.txt' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
PassedProbes_fa_out = '%(dir)s%(FirstG)sProbes_to_%(TarID)s_PASSED.fasta'%{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
MDfilled_sam_file = '%(dir)s%(FirstG)sProbes_aligned_to_%(TarID)s_MDfilled.sam'%{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
SamForTilingDepth = '%(dir)s%(FirstG)sProbes_PASSED_for%(TarID)s_sorted.sam' %{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
TilingDepth_out_filename = '%(dir)s%(FirstG)sProbes_PASSED_for%(TarID)s_sorted_tilingDepth.out' %{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
Gap_filled_coverage_out_filename = '%(dir)s%(FirstG)sProbes_PASSED_for%(TarID)s_sorted_tilingDepth_GapFiled.out'%{"dir":dirSubOuts,"FirstG":FirstGenomeID, "TarID":TargetID}
#HMR search for Target genome
OutPsl_Tar = '%(dir)sblat_%(TarID)s_to_human.psl' %{"dir":dirSubOuts, "TarID":TargetID}
HitsFasta_Tar_filename = '%(dir)sblat_%(TarID)s_to_human_BlockSizeFiltered_HitSeqs.fasta' %{"dir":dirSubOuts, "TarID":TargetID}
TargetSeqs_Tar = '%(dir)sblat_%(TarID)s_to_human_BlockSizeFiltered_HitSeqs.fasta'%{"dir":dirSubOuts, "TarID":TargetID}
OutBlast8_Tar_filename = '%(dir)sblat_%(TarID)s_to_human_BlockSizeFiltered.blast8' %{"dir":dirSubOuts, "TarID":TargetID}
filtered_out_Tar_filename = '%(dir)sblat_%(TarID)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_blast8.bed'%{"dir":dirSubOuts, "TarID":TargetID}
hmr_bedForTarGet_filename = '%(dir)sblat_%(TarID)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_HMRs.bed' %{"dir":dirSubOuts, "TarID":TargetID}
#New Probes out file in SureDesign format
New_probes_txt_out_filename = '%(dir)sNew_Probes_for_%(TarID)s.txt' %{"dir":dirSubOuts, "TarID":TargetID}
#New Probes out file in Fasta format
New_probes_fa_out_filename = '%(dir)sNew_Probes_for_%(TarID)s.fasta' %{"dir":dirSubOuts, "TarID":TargetID}
Updated_ProbeSetFasta = '%(dir)sUpdated_Probes.fasta' %{"dir":dirSubOuts}
Updated_ProbeSetTXT = '%(dir)sUpdated_Probes.txt' %{"dir":dirSubOuts}

#Step 2-1.
#build bwt index for G2
build_bwt_index(TargetGenomeFasta_name, Index_base)
#map ProbeSet against G2
bowtie_single(Index_base, TargetGenomeFasta_name, ProbesFasta, Out_sam_name)
#filter aligned Probes using mismatch, indel thresholds
Probe_aln_filter(ProbesTxtFile, MDfilled_sam_file, InDelThreshold, PassedProbes_fa_out, Out_sam_name, SamForTilingDepth, TotalIndelMiss, MaxPerMatch)
#Measure tiling depth
samtools_depth(TargetGenomeFasta_name, SamForTilingDepth, TilingDepth_out_filename, Gap_filled_coverage_out_filename)
#get HMR
######################cProfile.run('get_HMRs(TargetGenomeFasta_name, HitsFasta_Tar_filename, OutPsl_Tar, 0.0001, 40.0, OutBlast8_Tar_filename, filtered_out_Tar_filename, hmr_bedForTarGet_filename)')
#design new Probes for gaps
new_probes_design(TargetGenomeFasta_name, hmr_bedForTarGet_filename, Gap_filled_coverage_out_filename, New_probes_txt_out_filename, New_probes_fa_out_filename, TargetID, X)
#add new probes to ProbeSet
#Update probe set by including (adding) new produced probes
update_probeSet(ProbesFasta, New_probes_fa_out_filename, Updated_ProbeSetFasta)
update_probeSet(ProbesTxtFile, New_probes_txt_out_filename, Updated_ProbeSetTXT)
#return updated ProbeSet

probe_counts = '%(dir)sCounts_for_Updated_Probes.txt' %{"dir":dirSubOuts}

Targets = ['Mutu', 'Akata', 'GD1']
for item in Targets:
    Target_2ID = str(item)
    
    TargetGenomeFasta_name = '%(dir)s%(Tar_2ID)s_ebv.fa' %{"dir":dirInput, "Tar_2ID":Target_2ID}
    Index_base = '%(dir)s%(Tar_2ID)s_ebv' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    ProbesFasta = '%(dir)sUpdated_Probes.fasta' %{"dir":dirSubOuts}#'%(dir)s%(FirstG)s_Probes.fasta' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
    Out_sam_name = '%(dir)sUpdated_Probes_aligned_to_%(Tar_2ID)s.sam' %{"dir":dirSubOuts,"Tar_2ID":Target_2ID}
    ProbesTxtFile = '%(dir)sUpdated_Probes.txt' %{"dir":dirSubOuts}#'%(dir)s%(FirstG)s_Probes.txt' %{"dir":dirSubOuts, "FirstG":FirstGenomeID}
    PassedProbes_fa_out = '%(dir)sUpdated_Probes_%(Tar_2ID)s_PASSED.fasta'%{"dir":dirSubOuts,"Tar_2ID":Target_2ID}
    MDfilled_sam_file = '%(dir)sUpdated_Probes_aligned_to_%(Tar_2ID)s_MDfilled.sam'%{"dir":dirSubOuts,"Tar_2ID":Target_2ID}
    SamForTilingDepth = '%(dir)sUpdated_Probes_PASSED_for%(Tar_2ID)s_sorted.sam' %{"dir":dirSubOuts,"Tar_2ID":Target_2ID}
    TilingDepth_out_filename = '%(dir)sUpdated_Probes_PASSED_for%(Tar_2ID)s_sorted_tilingDepth.out' %{"dir":dirSubOuts,"Tar_2ID":Target_2ID}
    Gap_filled_coverage_out_filename = '%(dir)sUpdated_Probes_PASSED_for%(Tar_2ID)s_sorted_tilingDepth_GapFiled.out'%{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    #HMR search for Target genome
    OutPsl_Tar = '%(dir)sblat_%(Tar_2ID)s_to_human.psl' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    HitsFasta_Tar_filename = '%(dir)sblat_%(Tar_2ID)s_to_human_BlockSizeFiltered_HitSeqs.fasta' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    TargetSeqs_Tar = '%(dir)sblat_%(Tar_2ID)s_to_human_BlockSizeFiltered_HitSeqs.fasta'%{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    OutBlast8_Tar_filename = '%(dir)sblat_%(Tar_2ID)s_to_human_BlockSizeFiltered.blast8' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    filtered_out_Tar_filename = '%(dir)sblat_%(Tar_2ID)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_blast8.bed'%{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    hmr_bedForTarGet_filename = '%(dir)sblat_%(Tar_2ID)s_to_human_BlockSizeFiltered_eVal-ScoreFiltered_HMRs.bed' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    #New Probes out file in SureDesign format
    New_probes_txt_out_filename = '%(dir)sNew_Probes_for_%(Tar_2ID)s.txt' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    #New Probes out file in Fasta format
    New_probes_fa_out_filename = '%(dir)sNew_Probes_for_%(Tar_2ID)s.fasta' %{"dir":dirSubOuts, "Tar_2ID":Target_2ID}
    Updated_ProbeSetFasta = '%(dir)sUpdated_Probes.fasta' %{"dir":dirSubOuts}
    Updated_ProbeSetTXT = '%(dir)sUpdated_Probes.txt' %{"dir":dirSubOuts}




    #Step 2-2.
    #build bwt index for G2
    build_bwt_index(TargetGenomeFasta_name, Index_base)
    #map ProbeSet against G2
    bowtie_single(Index_base, TargetGenomeFasta_name, ProbesFasta, Out_sam_name)
    #filter aligned Probes using mismatch, indel thresholds
    Probe_aln_filter(ProbesTxtFile, MDfilled_sam_file, InDelThreshold, PassedProbes_fa_out, Out_sam_name, SamForTilingDepth, TotalIndelMiss, MaxPerMatch)
    #Measure tiling depth
    samtools_depth(TargetGenomeFasta_name, SamForTilingDepth, TilingDepth_out_filename, Gap_filled_coverage_out_filename)
    #get HMR
    #######################cProfile.run('get_HMRs(TargetGenomeFasta_name, HitsFasta_Tar_filename, OutPsl_Tar, 0.0001, 40.0, OutBlast8_Tar_filename, filtered_out_Tar_filename, hmr_bedForTarGet_filename)')
    #design new Probes for gaps
    new_probes_design(TargetGenomeFasta_name, hmr_bedForTarGet_filename, Gap_filled_coverage_out_filename, New_probes_txt_out_filename, New_probes_fa_out_filename, Target_2ID, X)
    #add new probes to ProbeSet
    #Update probe set by including (adding) new produced probes
    update_probeSet(ProbesFasta, New_probes_fa_out_filename, Updated_ProbeSetFasta)
    update_probeSet(ProbesTxtFile, New_probes_txt_out_filename, Updated_ProbeSetTXT)
    #return updated ProbeSet

#Probes_count(Updated_ProbeSetTXT, probe_counts)

Boosted_probes_txt_out_filename = '%(dir)s%(bst_type)s_Boosted_Probes.txt' %{"dir":dirSubOuts, "bst_type":Boost_Type}

boost(Updated_ProbeSetTXT, Boost_Type, Pool_factor, Boosted_probes_txt_out_filename)
