#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
File: virgena_config_generator.py
Created Date: February 29th 2020
Author: ZL Deng <dawnmsg(at)gmail.com>
---------------------------------------
Last Modified: 29th February 2020 5:07:17 pm
'''

import click
from os import path


@click.command()
@click.argument("r1", type=click.Path(exists=True))
@click.argument("r2", type=click.Path(exists=True))
@click.argument("ref", type=click.Path(exists=True))
@click.argument("outdir", type=str)
@click.option('-l', '--insert', type=int, default=800, help='Insrtion length')
@click.option('-t', '--threads', type=int, default=8, help='The number of threads')
def config_generator(r1, r2, ref, outdir, insert, threads):
    (r1, r2, ref, outdir) = list(map(path.abspath, [r1, r2, ref, outdir]))
    out_config = '''<config>
    <Data>
        <pathToReads1>{r1}</pathToReads1>
        <pathToReads2>{r2}</pathToReads2>
        <InsertionLength>{insert}</InsertionLength>
    </Data>
    <Reference>{ref}</Reference>
    <OutPath>{outdir}</OutPath>
    <ThreadNumber>{threads}</ThreadNumber>
	<BatchSize>10000</BatchSize>
    <ReferenceSelector>
		<Enabled>false</Enabled>
        <UseMajor>false</UseMajor>
        <ReferenceMSA></ReferenceMSA>
        <PathToUsearch>vsearch</PathToUsearch>
        <UclustIdentity>0.98</UclustIdentity>
        <MinReadLength>50</MinReadLength>
        <MinContigLength>1000</MinContigLength>
        <Delta>0.05</Delta>
		<MaxNongreedyComponentNumber>5</MaxNongreedyComponentNumber>
        <MapperToMSA>
			<K>7</K>
			<pValue>0.01</pValue>
			<IndelToleranceThreshold>1.5</IndelToleranceThreshold>
			<RandomModelParameters>
				<Order>4</Order>
				<ReadNum>10000</ReadNum>
				<Step>10</Step>
			</RandomModelParameters>
        </MapperToMSA>
        <Graph>
            <MinReadNumber>5</MinReadNumber>
            <VertexWeight>10</VertexWeight>
			<SimilarityThreshold>0.5</SimilarityThreshold>
			<Debug>false</Debug>
        </Graph>
        <Debug>false</Debug>
    </ReferenceSelector>
    <Mapper>
		<K>5</K>
        <pValue>0.01</pValue>
		<IndelToleranceThreshold>1.25</IndelToleranceThreshold>
		<RandomModelParameters>
			<Order>4</Order>
			<ReadNum>1000</ReadNum>
			<Step>10</Step>
		</RandomModelParameters>
        <Aligner>
            <Match>2</Match>
            <Mismatch>-3</Mismatch>
            <GapOpenPenalty>5</GapOpenPenalty>
            <GapExtensionPenalty>2</GapExtensionPenalty>
        </Aligner>
    </Mapper>
    <ConsensusBuilder>
        <IdentityThreshold>0.9</IdentityThreshold>
        <CoverageThreshold>0</CoverageThreshold>
        <MinIntersectionLength>10</MinIntersectionLength>
        <MinTerminationReadsNumber>1</MinTerminationReadsNumber>
		<Reassembler>
			<IdentityThreshold>0.9</IdentityThreshold>
			<MinTerminatingSequenceCoverage>0</MinTerminatingSequenceCoverage>
			<PairReadTerminationThreshold>0.1</PairReadTerminationThreshold>
			<MinReadLength>50</MinReadLength>
		</Reassembler>
		<Debug>false</Debug>
    </ConsensusBuilder>
    <Postprocessor>
        <Enabled>true</Enabled>
        <MinFragmentLength>500</MinFragmentLength>
        <MinIdentity>0.99</MinIdentity>
		<MinFragmentCoverage>0.99</MinFragmentCoverage>
        <Debug>false</Debug>
    </Postprocessor>
</config>'''.format(r1=r1, r2=r2, ref=ref, outdir=outdir, insert=insert, threads=threads)

    print(out_config)


if __name__ == '__main__':
    config_generator()
