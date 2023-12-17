#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========
serves as the wrapper for all guideseq pipeline

"""
import time
import os
import sys
import yaml
import argparse
import traceback

# Set up logger
import log
logger = log.createCustomLogger('root')

from alignReads import alignReads, index
from filterBackgroundSites import filterBackgroundSites
from umi import demultiplex, umitag, consolidate
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation

import multiprocessing

DEFAULT_DEMULTIPLEX_MIN_READS = 10000
DEFAULT_WINDOW_SIZE = 25
DEFAULT_MAX_SCORE = 7

CONSOLIDATE_MIN_QUAL = 15
CONSOLIDATE_MIN_FREQ = 0.9

def consolidate_util(fastq_file_read1, consolidated_fastq_file_read1, fastq_file_read2,
                         consolidated_fastq_file_read2, min_qual, min_freq):
        consolidate.consolidate(fastq_file_read1, consolidated_fastq_file_read1, min_qual, min_freq)
        consolidate.consolidate(fastq_file_read2, consolidated_fastq_file_read2, min_qual, min_freq)


class GuideSeq:
    def __init__(self):
        pass

    def parseManifest(self, manifest_path, without_contorl):
        logger.info('Loading manifest...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f)

        try:
            # Validate manifest data
            validation.validateManifest(manifest_data)

            self.BWA_path = manifest_data['bwa']
            self.bedtools = manifest_data['bedtools']
            self.reference_genome = manifest_data['reference_genome']
            self.output_folder = manifest_data['output_folder']
            self.undemultiplexed = manifest_data['undemultiplexed']
            self.samples = manifest_data['samples']

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS
        # Allow the user to specify window size for off-target search
        if 'window_size' in manifest_data:
            self.window_size = manifest_data['window_size']
        else:
            self.window_size = DEFAULT_WINDOW_SIZE
        # Allow the user to specify window size for off-target search
        if 'max_score' in manifest_data:
            self.max_score = manifest_data['max_score']
        else:
            self.max_score = DEFAULT_MAX_SCORE
        # Allow the user to specify PAM seq. Yichao 3/6/2020
        if 'PAM' in manifest_data:
            self.PAM = manifest_data['PAM']
        else:
            self.PAM = "NGG"

        if without_contorl:
            if len(self.samples) < 1:
                raise AssertionError('Your manifest does not contain any samples.')
        else:
            # Make sure the user has specified a control barcode
            if 'control' not in self.samples.keys():
                raise AssertionError('Your manifest must have a control sample specified.')

            # Make sure the user has both a sample and a control
            if len(self.samples) < 2:
                raise AssertionError('Your manifest must have at least one control and one treatment sample.')

        logger.info('Successfully loaded manifest.')

    def parseManifestDemultiplex(self, manifest_path):
        logger.info('Loading manifest for demultiplexing...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f)

            try:
                self.output_folder = manifest_data['output_folder']
                self.undemultiplexed = manifest_data['undemultiplexed']
                self.samples = manifest_data['samples']

            except Exception as e:
                logger.error('Incomplete or incorrect manifest file. Please ensure your manifest contains all required fields.')
                quit()

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS

        logger.info('Successfully loaded manifest for single-step demultiplexing.')

    def demultiplex(self):
        logger.info('Demultiplexing undemultiplexed files...')

        # Take our two barcodes and concatenate them
        swapped_sample_barcodes = {}
        for sample in self.samples:
            barcode1 = self.samples[sample]['barcode1']
            barcode2 = self.samples[sample]['barcode2']
            barcode = barcode1[1:8] + barcode2[1:8]
            swapped_sample_barcodes[barcode] = sample

        try:
            demultiplex.demultiplex(self.undemultiplexed['forward'],
                                    self.undemultiplexed['reverse'],
                                    self.undemultiplexed['index1'],
                                    self.undemultiplexed['index2'],
                                    swapped_sample_barcodes,
                                    os.path.join(self.output_folder, 'demultiplexed'),
                                    min_reads=self.demultiplex_min_reads)

            self.demultiplexed = {}
            for sample in self.samples:
                self.demultiplexed[sample] = {}
                self.demultiplexed[sample]['read1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r1.fastq')
                self.demultiplexed[sample]['read2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r2.fastq')
                self.demultiplexed[sample]['index1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i1.fastq')
                self.demultiplexed[sample]['index2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i2.fastq')

            logger.info('Successfully demultiplexed reads.')
        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            quit()

    def umitag(self, num_processes=20):
        logger.info('umitagging reads...')
        try:
            pool = multiprocessing.Pool(processes=num_processes)

            self.umitagged = {}
            for sample in self.samples:
                self.umitagged[sample] = {}
                self.umitagged[sample]['read1'] = os.path.join(self.output_folder, 'umitagged', sample + '.r1.umitagged.fastq')
                self.umitagged[sample]['read2'] = os.path.join(self.output_folder, 'umitagged', sample + '.r2.umitagged.fastq')

                job = pool.apply_async(umitag.umitag, args=(
                    self.demultiplexed[sample]['read1'],
                    self.demultiplexed[sample]['read2'],
                    self.demultiplexed[sample]['index1'],
                    self.demultiplexed[sample]['index2'],
                    self.umitagged[sample]['read1'],
                    self.umitagged[sample]['read2'],
                    os.path.join(self.output_folder, 'umitagged')))
                job.get()
                time.sleep(15) # Sleep for 15 seconds

            pool.close()
            pool.join()
            logger.info('Successfully umitagged reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def consolidate(self, min_freq=CONSOLIDATE_MIN_FREQ, min_qual=CONSOLIDATE_MIN_QUAL, num_processes=20):
        logger.info('Consolidating reads...')
        try:
            pool = multiprocessing.Pool(processes=num_processes)

            self.consolidated = {}

            for sample in self.samples:
                print(sample)
                self.consolidated[sample] = {}
                self.consolidated[sample]['read1'] = os.path.join(self.output_folder, 'consolidated', sample + '.r1.consolidated.fastq')
                self.consolidated[sample]['read2'] = os.path.join(self.output_folder, 'consolidated', sample + '.r2.consolidated.fastq')

                job = pool.apply_async(consolidate_util, args=(
                    self.umitagged[sample]['read1'], self.consolidated[sample]['read1'],
                    self.umitagged[sample]['read2'], self.consolidated[sample]['read2'],
                    min_qual, min_freq))
                job.get()
                time.sleep(15) # Sleep for 15 seconds

            pool.close()
            pool.join()
            logger.info('Successfully consolidated reads.')
        except Exception as e:
            logger.error('Error consolidating')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self, num_processes=20):
        logger.info('Aligning reads...')
        try:
            index(self.reference_genome, self.BWA_path)
            pool = multiprocessing.Pool(processes=num_processes)
            
            self.aligned = {}
            for sample in self.samples:
                sample_alignment_path = os.path.join(self.output_folder, 'aligned', sample + '.sam')
                self.aligned[sample] = sample_alignment_path
                job = pool.apply_async(alignReads, args=(
                    self.BWA_path, self.reference_genome,
                    self.consolidated[sample]['read1'],
                    self.consolidated[sample]['read2'],
                    sample_alignment_path))
                job.get()
                time.sleep(15) # Sleep for 15 seconds

            pool.close()
            pool.join()
            logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def identifyOfftargetSites(self):
        logger.info('Identifying offtarget sites...')

        try:
            self.identified = {}
            # Identify offtarget sites for each sample
            for sample in self.samples:

                # Prepare sample annotations
                sample_data = self.samples[sample]
                annotations = {}
                annotations['Description'] = sample_data['description']
                annotations['Targetsite'] = sample

                if sample is 'control':
                    annotations['Sequence'] = ''
                else:
                    annotations['Sequence'] = sample_data['target']

                samfile = os.path.join(self.output_folder, 'aligned', sample + '.sam')

                self.identified[sample] = os.path.join(self.output_folder, 'identified', sample + '_identifiedOfftargets.txt')

                identifyOfftargetSites.analyze(samfile, self.reference_genome, self.identified[sample], annotations,
                                               self.window_size, self.max_score)

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()

    def filterBackgroundSites(self):
        logger.info('Filtering background sites')

        try:
            self.filtered = {}

            # Filter background in each sample
            for sample in self.samples:
                if sample != 'control':
                    self.filtered[sample] = os.path.join(self.output_folder, 'filtered', sample + '_backgroundFiltered.txt')
                    filterBackgroundSites(self.bedtools, self.identified[sample], self.identified['control'], self.filtered[sample])
                    logger.info('Finished background filtering for {0} sample'.format(sample))

            logger.info('Finished filtering background sites.')

        except Exception as e:
            logger.error('Error filtering background sites.')
            logger.error(traceback.format_exc())

    def visualize(self):
        logger.info('Visualizing off-target sites')
        for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
            if sample != 'control':
                try:
                    infile = self.identified[sample]
                    outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')
                    try:
                        self.PAM
                        visualizeOfftargets(infile, outfile, title=sample,PAM=self.PAM)
                    except:
                        visualizeOfftargets(infile, outfile, title=sample,PAM="NGG")
                except Exception as e:
                    logger.error('Error visualizing off-target sites: %s'%(sample))
                    logger.error(traceback.format_exc())
        logger.info('Finished visualizing off-target sites')


def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    all_parser.add_argument('--identifyAndFilter', action='store_true', default=False)
    all_parser.add_argument('--skip_demultiplex', action='store_true', default=False)
    all_parser.add_argument('--skip_umitag', action='store_true', default=False)
    all_parser.add_argument('--skip_consolidate', action='store_true', default=False)
    all_parser.add_argument('--skip_align', action='store_true', default=False)
    all_parser.add_argument('--without_contorl', action='store_true', default=False)
    all_parser.add_argument('--n_workers', default=1, type=int)

    demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplex undemultiplexed FASTQ files')
    demultiplex_parser.add_argument('--manifest', '-m', help='Specify the manifest path', required=True)

    umitag_parser = subparsers.add_parser('umitag', help='UMI tag demultiplexed FASTQ files for consolidation')
    umitag_parser.add_argument('--read1', required=True)
    umitag_parser.add_argument('--read2', required=True)
    umitag_parser.add_argument('--index1', required=True)
    umitag_parser.add_argument('--index2', required=True)
    umitag_parser.add_argument('--outfolder', required=True)

    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate UMI tagged FASTQs')
    consolidate_parser.add_argument('--read1', required=True)
    consolidate_parser.add_argument('--read2', required=True)
    consolidate_parser.add_argument('--outfolder', required=True)
    consolidate_parser.add_argument('--min_quality', default=CONSOLIDATE_MIN_QUAL, type=float)
    consolidate_parser.add_argument('--min_frequency', default=CONSOLIDATE_MIN_FREQ, type=float)
    consolidate_parser.add_argument('--n_workers', default=1, type=int)

    align_parser = subparsers.add_parser('align', help='Paired end read mapping to genome')
    align_parser.add_argument('--bwa', required=True)
    align_parser.add_argument('--genome', required=True)
    align_parser.add_argument('--read1', required=True)
    align_parser.add_argument('--read2', required=True)
    align_parser.add_argument('--outfolder', required=True)
    align_parser.add_argument('--n_workers', default=1, type=int)

    identify_parser = subparsers.add_parser('identify', help='Identify GUIDE-seq offtargets')
    identify_parser.add_argument('--aligned', required=True)
    identify_parser.add_argument('--genome', required=True)
    identify_parser.add_argument('--outfolder', required=True)
    identify_parser.add_argument('--target_sequence', required=True)
    identify_parser.add_argument('--description', required=False)
    identify_parser.add_argument('--max_score', required=False, type=int, default=7)
    identify_parser.add_argument('--window_size', required=False, type=int, default=25)

    filter_parser = subparsers.add_parser('filter', help='Filter identified sites from control sites')
    filter_parser.add_argument('--bedtools', required=True)
    filter_parser.add_argument('--identified', required=True)
    filter_parser.add_argument('--background', required=True)
    filter_parser.add_argument('--outfolder', required=True)

    visualize_parser = subparsers.add_parser('visualize', help='Visualize off-target sites')
    visualize_parser.add_argument('--infile', required=True)
    visualize_parser.add_argument('--outfolder', required=True)
    visualize_parser.add_argument('--title', required=False)

    return parser.parse_args()


def load_demultiplexed(args):
    g = GuideSeq()
    g.parseManifest(args.manifest, args.without_contorl)
    g.demultiplexed = {}
    for sample in g.samples:
        g.demultiplexed[sample] = {}
        g.demultiplexed[sample]['read1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r1.fastq')
        g.demultiplexed[sample]['read2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r2.fastq')
        g.demultiplexed[sample]['index1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i1.fastq')
        g.demultiplexed[sample]['index2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i2.fastq')
        if not os.path.isfile(g.demultiplexed[sample]['read1']):
            print ("Can't find ", g.demultiplexed[sample]['read1'])
            exit()
        if not os.path.isfile(g.demultiplexed[sample]['read2']):
            print ("Can't find ", g.demultiplexed[sample]['read2'])
            exit()
        if not os.path.isfile(g.demultiplexed[sample]['index1']):
            print ("Can't find ", g.demultiplexed[sample]['index1'])
            exit()
        if not os.path.isfile(g.demultiplexed[sample]['index2']):
            print ("Can't find ", g.demultiplexed[sample]['index2'])
            exit()

    return g


def load_umitagged(g):
    g.umitagged = {}
    for sample in g.samples:
        g.umitagged[sample] = {}
        g.umitagged[sample]['read1'] = os.path.join(g.output_folder, 'umitagged', sample + '.r1.umitagged.fastq')
        g.umitagged[sample]['read2'] = os.path.join(g.output_folder, 'umitagged', sample + '.r2.umitagged.fastq')
        if not os.path.isfile(g.umitagged[sample]['read1']):
            print ("Can't find ", g.umitagged[sample]['read1'])
            exit()
        if not os.path.isfile(g.umitagged[sample]['read2']):
            print ("Can't find ", g.umitagged[sample]['read2'])
            exit()


def load_consolidate(g):
    g.consolidated = {}
    for sample in g.samples:
        print(sample)
        g.consolidated[sample] = {}
        g.consolidated[sample]['read1'] = os.path.join(g.output_folder, 'consolidated', sample + '.r1.consolidated.fastq')
        g.consolidated[sample]['read2'] = os.path.join(g.output_folder, 'consolidated', sample + '.r2.consolidated.fastq')
        if not os.path.isfile(g.consolidated[sample]['read1']):
            print ("Can't find ", g.consolidated[sample]['read1'])
            exit()
        if not os.path.isfile(g.consolidated[sample]['read2']):
            print ("Can't find ", g.consolidated[sample]['read2'])
            exit()

def load_align(g):
    g.aligned = {}
    for sample in g.samples:
        sample_alignment_path = os.path.join(g.output_folder, 'aligned', sample + '.sam')
        g.aligned[sample] = sample_alignment_path
        if not os.path.isfile(g.aligned[sample]):
            print ("Can't find ", g.aligned[sample])
            exit()


def main():
    args = parse_args()
    # print("waiting some time")
    # time.sleep(3600*10)
    # print("end wait time")
    if args.command == 'all':
        if args.identifyAndFilter:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest, args.without_contorl)
                load_align(g)

                g.identifyOfftargetSites()
                if args.without_contorl:
                     logger.info('skipping filter background sites as not control is provided.')
                else:
                    g.filterBackgroundSites()
                g.visualize()
            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_demultiplex:
            try:
                g = load_demultiplexed(args)
                g.umitag(num_processes=args.n_workers)
                g.consolidate(num_processes=args.n_workers)
                g.alignReads(num_processes=args.n_workers)
                g.identifyOfftargetSites()
                if args.without_contorl:
                     logger.info('skipping filter background sites as not control is provided.')
                else:
                    g.filterBackgroundSites()
                g.visualize()
            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_umitag:
            try:
                g = load_demultiplexed(args)
                load_umitagged(g)
                g.consolidate(num_processes=args.n_workers)
                g.alignReads(num_processes=args.n_workers)
                g.identifyOfftargetSites()
                if args.without_contorl:
                     logger.info('skipping filter background sites as not control is provided.')
                else:
                    g.filterBackgroundSites()
                g.visualize()
            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_consolidate:
            try:
                g = load_demultiplexed(args)
                load_umitagged(g)
                load_consolidate(g)
                g.alignReads(num_processes=args.n_workers)
                g.identifyOfftargetSites()
                if args.without_contorl:
                    logger.info('skipping filter background sites as not control is provided.')
                else:
                    g.filterBackgroundSites()
                g.visualize()
            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_align:
            try:
                g = load_demultiplexed(args)
                load_umitagged(g)
                load_consolidate(g)
                load_align(g)
                g.identifyOfftargetSites()
                if args.without_contorl:
                    logger.info('skipping filter background sites as not control is provided.')
                else:
                    g.filterBackgroundSites()
                g.visualize()
            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        else:
            g = GuideSeq()
            g.parseManifest(args.manifest, args.without_contorl)
            g.demultiplex()
            g.umitag(num_processes=args.n_workers)
            g.consolidate(num_processes=args.n_workers)
            g.alignReads(num_processes=args.n_workers)
            g.identifyOfftargetSites()
            if args.without_contorl:
                logger.info('skipping filter background sites as not control is provided.')
            else:
                g.filterBackgroundSites()
            g.visualize()
    elif args.command == 'demultiplex':
        """
        Run just the demultiplex step given the manifest
        """
        g = GuideSeq()
        g.parseManifestDemultiplex(args.manifest)
        g.demultiplex()
    elif args.command == 'umitag':
        """
        Run just the umitag step
        python guideseq/guideseq.py umitag --read1 test/data/demultiplexed/EMX1.r1.fastq --read2 test/data/demultiplexed/EMX1.r2.fastq --index1 test/data/demultiplexed/EMX1.i1.fastq --index2 test/data/demultiplexed/EMX1.i2.fastq --outfolder test/output/
        """
        g = GuideSeq()
        g.output_folder = args.outfolder
        sample = os.path.basename(args.read1).split('.')[0]
        g.samples = [sample]
        g.demultiplexed = {sample: {}}
        g.demultiplexed[sample]['read1'] = args.read1
        g.demultiplexed[sample]['read2'] = args.read2
        g.demultiplexed[sample]['index1'] = args.index1
        g.demultiplexed[sample]['index2'] = args.index2
        g.umitag(num_processes=args.n_workers)
    elif args.command == 'consolidate':
        """
        Run just the consolidate step
        python guideseq/guideseq.py consolidate --read1 test/data/umitagged/EMX1.r1.umitagged.fastq --read2 test/data/umitagged/EMX1.r2.umitagged.fastq --outfolder test/output/ --min_frequency 0.8 --min_quality 14
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.umitagged = {sample: {}}
        g.umitagged[sample]['read1'] = args.read1
        g.umitagged[sample]['read2'] = args.read2
        g.consolidate(min_freq=args.min_frequency, min_qual=args.min_quality, num_processes=args.n_workers)
    elif args.command == 'align':
        """
        Run just the alignment step
        python guideseq/guideseq.py align --bwa bwa --read1 test/data/consolidated/EMX1.r1.consolidated.fastq --read2 test/data/consolidated/EMX1.r2.consolidated.fastq --genome /Volumes/Media/hg38/hg38.fa --outfolder test/output/
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.BWA_path = args.bwa
        g.reference_genome = args.genome
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.consolidated = {sample: {}}
        g.consolidated[sample]['read1'] = args.read1
        g.consolidated[sample]['read2'] = args.read2
        g.alignReads(num_processes=args.n_workers)
    elif args.command == 'identify':
        """
        Run just the identify step
        python guideseq/guideseq.py identify --genome /Volumes/Media/hg38/hg38.fa --aligned test/output/aligned/EMX1.sam --outfolder test/output/ --target_sequence GAGTCCGAGCAGAAGAAGAANGG
        """
        if 'description' in args:
            description = args.description
        else:
            description = ''

        if 'max_score' in args:
            max_score = args.max_score
        else:
            max_score = 7

        if 'window_size' in args:
            window_size = args.window_size
        else:
            window_size = 25

        g = GuideSeq()
        g.output_folder = args.outfolder
        g.reference_genome = args.genome
        sample = os.path.basename(args.aligned).split('.')[0]
        g.samples = {sample: {'description': description, 'target': args.target_sequence}}
        g.aligned = {sample: args.aligned}
        g.max_score = max_score
        g.window_size = window_size
        g.identifyOfftargetSites()

    elif args.command == 'filter':
        """
        Run just the filter step

        """
        sample = os.path.basename(args.identified).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.bedtools = args.bedtools
        g.samples = {sample: {}, 'control': {}}
        g.identified = {}
        g.identified[sample] = args.identified
        g.identified['control'] = args.background
        g.filterBackgroundSites()
    elif args.command == 'visualize':
        """
        Run just the visualize step
        """
        g = GuideSeq()
        g.output_folder = os.path.dirname(args.outfolder)
        sample = os.path.basename(args.infile).split('.')[0]
        g.samples = {sample: {}}
        g.identified = {}
        g.identified[sample] = args.infile
        g.visualize()


if __name__ == '__main__':
    main()
