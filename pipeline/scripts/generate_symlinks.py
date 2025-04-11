#!/usr/bin/env python

import subprocess
import argparse

import os
import shutil
import logging
import re
from collections import defaultdict


class GenerateSymlink():
    "Generates symlinks required for IGVnav inputs"
    
    def __init__(self, outputdirname, scriptdirname, targets):
        """ 
        :params: outputdirname -  sample output directory name 
        """
        self.outputdirname = outputdirname
        self.script_dir = scriptdirname
        self.targets = targets
        self.is_wgs = True if targets == '' else False
    

    def generateIGVsymlink(self, *args):
        """
        Generates the symlinks for give list of tuples
        :params: args - ('directoryname', 'patters')
        """   
        igvnav_dirname_dst = os.path.join(self.outputdirname, 'IGVnav')
        src_dir = os.path.abspath(self.outputdirname)
        logging.info("Generating IGVNav Symlinks in : " + igvnav_dirname_dst ) 
        if not os.path.exists(igvnav_dirname_dst): os.mkdir(igvnav_dirname_dst)
        try:
            if self.is_wgs:
                symlinks = (('variants','.vep.vcf.gz'),('bams','nodups.bam'), ('bams','nodups.bam.bai'),
                            ('bams', 'markdups.bam'), ('bams', 'markdups.bam.bai'),
                            ('bams','clipoverlap.bam'), ('bams','clipoverlap.bai'), ('cnv', '.bedGraph'),
                            ('variants', '.bedGraph'), ('svs/igv','.mut'), ('svs','.gtf'), 
                            ('svs/gridss', 'evidence.bam'), ('svs/gridss', 'evidence.bam.bai'), 
                            ('svs/gridss', 'assembly.bam'), ('svs','.bam'), ('svs','.bai'),
                            ('', 'igvnav-input.txt')) + args
            else:
                symlinks = (('variants','.vep.vcf'),('bams','nodups.bam'), ('bams','nodups.bam.bai'),
                            ('bams','clipoverlap.bam'), ('bams','clipoverlap.bai'), ('variants','.vep.vcf'), 
                            ('cnv', '.bedGraph'), ('variants', '.bedGraph'), ('svs/igv','.mut'), ('svs','.gtf'), 
                            ('svs/gridss', 'evidence.bam'), ('svs/gridss', 'evidence.bam.bai'), 
                            ('svs','.bam'), ('svs','.bai'), ('', 'igvnav-input.txt')) + args
            
            for each_input in symlinks:
                dir_name = os.path.join(src_dir,each_input[0])
                self.create_symlink(dir_name, src_dir, igvnav_dirname_dst, each_input[1])
            logging.info("Created IGVNav Symlinks : " + igvnav_dirname_dst ) 
        except Exception as e:
            logging.info(e)
        return  

    def create_symlink(self, travers_dir_name, src_dir, igvnav_dirname_dst, suffix):
        "Recursively Traverse through the directory and create symlink"
        for root, dirs, files in os.walk(travers_dir_name):
            for each_file in files:
                if each_file.endswith(suffix) and not os.path.exists(os.path.join(igvnav_dirname_dst,each_file)):
                    try:		
                        os.symlink(os.path.join(root,each_file), os.path.join(igvnav_dirname_dst,each_file))
                    except Exception as e:
                        logging.info(e) 
        return

    def get_all_files(self, dir_path):
        igv_session_file_list = defaultdict(lambda : defaultdict(list))
        all_files = [('bam_common', 'bam_nodups', '.*nodups.bam$'),
                     ('sv', 'mut_svaba_somatic','.*_(somatic)_svaba.mut$'),
                     ('sv', 'mut_gridss_cfdna', '.*-CFDNA-.*_gridss.mut$'),
                     ('sv', 'mut_svcaller_cfdna', '.*-CFDNA-.*_svcaller.mut$'),
                     ('sv', 'mut_lumpy', '.*lumpy_len500_SU24.mut$'),
                     ('sv', 'mut_svaba_germline','.*_(germline)_svaba.mut$'),
                     ('sv', 'mut_gridss_normal', '^(?:(?!CFDNA).)*_gridss.mut$'),
                     ('sv', 'mut_svcaller_normal', '^(?:(?!CFDNA).)*_svcaller.mut$'),
                     ('sv', 'gtf_cfdna', '.*-CFDNA-.*svs.gtf$'),
                     ('sv', 'gtf_normal', '^(?:(?!CFDNA).)*svs.gtf$'),
                     ('snps', 'bam_cfdna', '.*-CFDNA-.*_clipoverlap.bam$'),
                     ('snps', 'bam_normal', '^(?:(?!CFDNA).)*clipoverlap.bam$'),
                     ('cnv', 'flank_profile_cfdna', '.*-CFDNA-.*_profile.bedGraph'),
                     ('cnv', 'flank_profile_normal', '^(?:(?!CFDNA).)*_profile.bedGraph'),
                     ('cnv', 'flank_cnv_cfdna', '.*-CFDNA-.*_segments.bedGraph'),
                     ('cnv', 'flank_cnv_normal', '^(?:(?!CFDNA).)*_segments.bedGraph'),
                     ('asf', 'flank_asf', '.*germline-variants-somatic-afs.bedGraph'),
                     ]
            
        if self.is_wgs:
            all_files.extend([('bam_hmftools', 'bam_markdups', '.*markdups.bam$'),
                              ('snps', 'vep', '.*.all.(somatic|germline).vep.vcf.gz$'),
                              ('sv', 'bam_gridss_normal', '^(?:(?!CFDNA|T).)*evidence.bam$'),
                              ('sv', 'bam_gridss_tumor', '.*-(T|CFDNA).*evidence.bam$')
                            ])
        else:
            all_files.extend([('sv', 'bam_cfdna', '.*-CFDNA-.*(svs|evidence).bam$'),
                              ('sv', 'bam_normal', '^(?:(?!CFDNA).)*(svs|evidence).bam$'),
                              ('snps', 'vep', '.*.all.(somatic|germline).vep.vcf$')
                            ])


        for variant_type, file_type, regex in all_files:
            files_list = list(filter( lambda x: re.match(regex, x) and not x.startswith('.') and not x.endswith('.out'), os.listdir(dir_path)))
            files_list.sort()
            igv_session_file_list[variant_type][file_type]=files_list

        return igv_session_file_list

    def create_igv_session_file(self):
        """
        Create IGV seesion file in xml format for given sample
        """

        igv_session_master=os.path.join(self.script_dir, "igv_session_master.xml")
        igv_session_sv_master=os.path.join(self.script_dir, "igv_session_sv_master.xml")

        resource_path = """
        <Resource path="{path_name}"/>
        """
        panel_str = """
        <Panel height="{panel_height}"  width="{panel_width}" name="{panel_name}">
            {tracks}
        </Panel>
        """
        tracks_bam_str= """
        <Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="{color}" colorScale="ContinuousColorScale;0.0;927.0;255,255,255;{color}" fontSize="10" id="{bam_file_full}_coverage" name="{bam_file} Coverage" snpThreshold="0.0001" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1503.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" height="60" id="{bam_file_full}_junctions" name="{bam_file} Junctions" visible="false"/>
        <Track clazz="org.broad.igv.sam.AlignmentTrack" displayMode="EXPANDED" experimentType="OTHER" fontSize="10" id="{bam_file_full}" name="{bam_file}" visible="true">
            <RenderOptions colorOption="READ_STRAND"/>
        </Track>
        """
        
        track_capture_bed = """
        <Track clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;103.0;255,255,255;0,0,178" fontSize="10" id="{capture_bed_full}" name="{capture_bed}" visible="true"/>
        """
        
        sv_gtf_track_str="""
            <Track clazz="org.broad.igv.track.FeatureTrack" color="{color}" fontSize="10" id="{gtf_file_full_path}" name="{gtf_file_path}" visible="true"/>
        """
        sv_mut_track_str="""
                <Track clazz="org.broad.igv.track.MutationTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;12.0;255,255,255;0,0,178" fontSize="10" height="15" id="{mut_file_full_path}_{id_field}" name="{id_field}" visible="true"/>
        """
        snp_vcf = """
            <Track clazz="org.broad.igv.variant.VariantTrack" color="0,0,178" displayMode="EXPANDED" fontSize="10" id="{vcf_full_file_path}" name="{vcf_file_path}" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>
        """
        cnv_flank = """
            <Track autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{full_path}" name="{file_name}" altColor="0,0,255" color="255,0,0"  renderer="BAR_CHART" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.0" minimum="-1.0" type="LINEAR"/>
            </Track>
        """
        cnv_flank_profile = """
            <Track  autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{full_path}" name="{file_name}" altColor="0,0,255" color="255,0,0" renderer="SCATTER_PLOT" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.0" minimum="-1.0" type="LINEAR"/>
        </Track>

        """
        snp_asf = """
            <Track  clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{full_path}" name="SNP_tDNA"  renderer="SCATTER_PLOT" color="0,0,0"  visible="true" windowFunction="mean">
            <DataRange baseline="0.5" drawBaseline="true" flipAxis="false" maximum="1" minimum="0" type="LINEAR"/>
        </Track>

        """


        # try:
        logging.info(" Generating IGV Session File ")
        igvnav_dir = os.path.join(self.outputdirname, 'IGVnav')
        if not os.path.exists(igvnav_dir): raise Exception('IGVnav folder not found..')
        if not os.path.exists(igv_session_master): raise Exception('IGV session master file not found : /nfs/PROBIO/for_igv/igv_session_master.xml ' )
        if not os.path.exists(igv_session_master): raise Exception('IGV session master file not found : /nfs/PROBIO/for_igv/igv_session_sv_master.xml ' )
        igv_session_master_str=open(igv_session_master, 'r').read()
        igv_session_sv_master_str=open(igv_session_sv_master, 'r').read()
        igv_session_files = self.get_all_files(igvnav_dir)
        
        # target_name = capture_kit_loopkup[capture_id[0:2]]
        target_bed = self.targets

        #session file for snps
        snp_resource = ""
        snp_bam_panel = ""
        snp_vep = ""
        type_color_arr = {'DEL':'53,116,199', 'TRA': '239,133,62', 'DUP' :'204,55,48', 'INV' :'69,136,51', 'COM' : '2,2,0'}

        capture_bed = ''
        if target_bed != '':
            snp_resource += resource_path.format(path_name=target_bed)
            capture_bed += track_capture_bed.format(capture_bed_full=target_bed, capture_bed=os.path.basename(target_bed))

        all_snp_files = [('snps', 'bam_cfdna'), ('snps', 'bam_normal'),
                         ('cnv', 'flank_cnv_normal'),
                         ('cnv', 'flank_profile_normal'),
                         ('cnv', 'flank_cnv_cfdna'),
                         ('cnv', 'flank_profile_cfdna'),
                         ('asf', 'flank_asf'),
                         ('snps', 'vep')]

        if self.is_wgs:
            all_snp_files.extend([('bam_common', 'bam_nodups'),
                                  ('bam_hmftools', 'bam_markdups')])

        for track_type, each_track in all_snp_files:
            for each_file in igv_session_files[track_type][each_track]:
                full_path = igvnav_dir + '/' + each_file
                if not os.path.exists(full_path) or not os.path.getsize(full_path):
                    continue
                ### Skip vep.vcf.gz file for WGS because it is taking longer time to load 
                ### Due to millions of records
                ### TODO filter out variants from vep.vcf.gz file
                if self.is_wgs and each_track.startswith('vep'):
                    continue
                snp_resource += resource_path.format(path_name=full_path)
                if each_track.startswith('bam'):
                    regex = ".*-(N|CFDNA|T)-.*(DEL|DUP|INV|TRA).bam$"
                    matches = re.search(regex, each_file)
                    if matches:
                        bam_color = str(type_color_arr[matches.group(2)])
                    else:
                        bam_color = '175,175,175'
                    snp_bam_track = tracks_bam_str.format(bam_file_full=full_path, bam_file=each_file, color=bam_color)
                    snp_bam_panel += panel_str.format(panel_height=250, panel_width=1800, panel_name=each_file , tracks=snp_bam_track)
                if each_track.startswith('flank_cnv'):
                    file_path_arr = full_path.split('/')[-1]
                    if  re.match('^(?:(?!-(CFDNA|T)-).)*_segments.bedGraph', file_path_arr):
                        file_name = 'gDNA_segments'
                    elif re.match('.*-(CFDNA|T)-.*_segments.bedGraph', file_path_arr):
                        file_name = 'tDNA_segments'
                    else:
                        file_name = full_path
                    snp_vep += cnv_flank.format(full_path=full_path, file_name=file_name)
                if each_track.startswith('flank_profile'):
                    file_path_arr = full_path.split('/')[-1]
                    if re.match('^(?:(?!-(CFDNA|T)-).)*_profile.bedGraph', file_path_arr):
                        file_name = 'gDNA_bins'
                    elif re.match('.*-(CFDNA|T)-.*_profile.bedGraph', file_path_arr):
                        file_name = 'tDNA_bins'
                    else:
                        file_name = full_path
                    snp_vep += cnv_flank_profile.format(full_path=full_path, file_name=file_name)
                if each_track.startswith('flank_asf'):
                    snp_vep += snp_asf.format(full_path=full_path)
                if each_track.startswith('vep'):
                    file_split = each_file.split('-all.')
                    vcf_file_name = file_split[1].replace('.',' ').upper()
                    snp_vep += snp_vcf.format(vcf_full_file_path=full_path, vcf_file_path=vcf_file_name)

        #session file for structural variants
        all_sv_files=[('bam_common', 'bam_nodups'), ('bam_hmftools', 'bam_markdups'),
                      ('sv', 'bam_cfdna'), ('sv', 'bam_normal'), 
                      ('cnv', 'flank_cnv_normal'), ('cnv', 'flank_profile_normal'),
                      ('cnv', 'flank_cnv_cfdna'), ('cnv', 'flank_profile_cfdna'),
                      ('asf', 'flank_asf'), ('sv', 'mut_svaba_somatic'), 
                      ('sv', 'mut_svcaller_cfdna'), ('sv', 'mut_gridss_cfdna'),
                      ('sv', 'mut_lumpy'), ('sv', 'gtf_cfdna'),
                      ('sv', 'mut_svaba_germline'), ('sv', 'mut_svcaller_normal'), 
                      ('sv', 'mut_gridss_normal'), ('sv', 'gtf_normal'), 
                      ('sv', 'bam_gridss_normal'), ('sv', 'bam_gridss_tumor')]

        sv_resource= ""
        sv_bam_panel= ""
        sv_mut_track= ""
        sv_gtf_track= ""

        if target_bed != '':
            sv_resource += resource_path.format(path_name=target_bed)

        for track_type, each_track in all_sv_files:
            for each_file in igv_session_files[track_type][each_track]:

                full_path = igvnav_dir + '/' + each_file
                if not os.path.exists(full_path) :
                    continue

                if not os.path.getsize(full_path):
                    continue

                #sv_resource += resource_path.format(path_name=full_path)
                
                if each_track.startswith('bam'):
                    sv_resource += resource_path.format(path_name=full_path)
                    regex = ".*-(N|CFDNA|T)-.*(DEL|DUP|INV|TRA).bam$"
                    matches = re.search(regex, each_file)
                    if matches:
                        bam_color = str(type_color_arr[matches.group(2)])
                    else:
                        bam_color = '175,175,175'
                    sv_bam_track = tracks_bam_str.format(bam_file_full=full_path, bam_file=each_file, color=bam_color)
                    sv_bam_panel += panel_str.format(panel_height=50, panel_width=2543, panel_name=each_file, tracks=sv_bam_track)

                if each_track.startswith('mut'):
                    # get content of SDID column to append to id and name fields in the track
                    f = open(full_path, 'r')
                    line = f.readline() # skip the header line
                    line = f.readline()  # get the first variant line
                    f.close()
                    if line:
                        id_field = line.strip().split()[3]  # use the 4th column (SDID column) if there is a variant
                    else:
                        continue  # skip the file if there are no variants

                    sv_resource += resource_path.format(path_name=full_path)
                    sv_mut_track += sv_mut_track_str.format(mut_file_full_path=full_path, id_field=id_field)

                if each_track.startswith('flank_cnv'):
                    file_path_arr = full_path.split('/')[-1]
                    if  re.match('^(?:(?!-(CFDNA|T)-).)*_segments.bedGraph', file_path_arr):
                        file_name = 'gDNA_segments'
                    elif re.match('.*-(CFDNA|T)-.*_segments.bedGraph', file_path_arr):
                        file_name = 'tDNA_segments'
                    else:
                        file_name = full_path
                    sv_resource += resource_path.format(path_name=full_path)
                    sv_mut_track += cnv_flank.format(full_path=full_path, file_name=file_name)
                if each_track.startswith('flank_profile'):
                    file_path_arr = full_path.split('/')[-1]
                    if re.match('^(?:(?!-(CFDNA|T)-).)*_profile.bedGraph', file_path_arr):
                        file_name = 'gDNA_bins'
                    elif re.match('.*-(CFDNA|T)-.*_profile.bedGraph', file_path_arr):
                        file_name = 'tDNA_bins'
                    else:
                        file_name = full_path
                    sv_resource += resource_path.format(path_name=full_path)
                    sv_mut_track += cnv_flank_profile.format(full_path=full_path, file_name=file_name)
                if each_track.startswith('flank_asf'):
                    sv_resource += resource_path.format(path_name=full_path)
                    sv_mut_track += snp_asf.format(full_path=full_path)

                if each_track.startswith('gtf'):
                    sv_resource += resource_path.format(path_name=full_path)
                    regex = ".*-(N|CFDNA|T)-.*(DEL|DUP|INV|TRA).gtf$"
                    matches = re.search(regex, each_file)
                    if matches:
                        gtf_color = str(type_color_arr[matches.group(2)])
                    else:
                        gtf_color = '0,0,178'

                    sv_gtf_track += sv_gtf_track_str.format(gtf_file_full_path=full_path, gtf_file_path=each_file, color=gtf_color)


        sv_session_data=igv_session_sv_master_str.format(add_resource=sv_resource, 
                                                            add_sv_mut_track=sv_mut_track, 
                                                            add_panel=sv_bam_panel,
                                                            add_capture_bed=capture_bed,
                                                            add_sv_gtf_track=sv_gtf_track)
        snp_session_data=igv_session_master_str.format(add_resource=snp_resource, 
                                                        add_panel=snp_bam_panel, 
                                                        add_capture_bed=capture_bed, 
                                                        add_vcf_track=snp_vep)
        
        with open(igvnav_dir+'/igv_session_snps.xml', 'w') as fw:
            fw.write(snp_session_data)

        with open(igvnav_dir+'/igv_session_sv.xml', 'w') as fw:
            fw.write(sv_session_data)

        logging.info("Created IGV Session Files..")

        # except Exception as e:
        #     print('error', e)
        #     logging.info("Error While creating session file : " + str(e))


if __name__ == "__main__":
    """
    create symlinks for IGVNav
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--targets', help="Target bed file - capture kit id")
    parser.add_argument('--outdir', required=True, help="project output dir")
    parser.add_argument('--script-dir', help="script dir for IGV session xml files")
    args = parser.parse_args()
    targets = ''
    if args.targets:
        targets = args.targets
    
    if not args.script_dir:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    else:
        script_dir = args.script_dir
    
    create_symlinks = GenerateSymlink(args.outdir, script_dir, targets)
    create_symlinks.generateIGVsymlink()
    create_symlinks.create_igv_session_file()
