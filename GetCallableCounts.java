/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//package pancancer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * modify Burair's Algorithms
 * @author Burair
 */
public class GetCallableCounts {
    
    public static void main(String[] args) throws IOException {
       
        String annotationFile = "G:\\MAC_Research_Data\\Pan_cancer\\Mapping_source\\Canis_familiaris.CanFam3.1.99.chr.gtf";
        String geneListFile = "G:\\MAC_Research_Data\\Pan_cancer\\Retro_gene_finding\\RetroGeneList\\new_retro_gene_list_CanFam3.1.99gtf.txt";
        String sampleWigFile = "G:\\MAC_Research_Data\\Pan_cancer\\Retro_gene_finding\\sample_wig.txt";
        String outputFile = "G:\\MAC_Research_Data\\Pan_cancer\\Retro_gene_finding\\new_GTF-java_test.txt";
        
        // Step 1: Reading gene list
        System.out.println("Step 1: reading gene list file: " + geneListFile);
        ArrayList<String> geneList = new ArrayList<>();
        Set<String> uniqueGenes = new TreeSet<>();
        BufferedReader reader = new BufferedReader(new FileReader(geneListFile));
        String line = reader.readLine();
        while(line != null) {
            if(line.trim().isEmpty() == false && line.startsWith("#") == false) {
                String gene = line.trim();
                if(uniqueGenes.add(gene) == true) {
                    geneList.add(gene);
                }
            }
            line = reader.readLine();
        }
        uniqueGenes = null; // not needed anymore
        reader.close();
        System.out.println("Finished reading gene list file: " + geneListFile);
        System.out.println(geneList.size() + " unique genes found");
        System.out.println();
        
        // Step 2: Reading annotation file (gene info)
        System.out.println("Step 2: reading annotation file: " + annotationFile);
        Map<String, GeneLocationInfo> geneInfo = getGenesInfo(annotationFile, geneList);
        System.out.println("Finished reading annotation file: " + annotationFile);
        System.out.println(geneInfo.size() + " of " + geneList.size() + " genes are annotated");
        System.out.println();
        
        // Step 3: Processing sample wig files
        System.out.println("Step 3: processing sample wig files ...");
        ArrayList<String> sampleNames = new ArrayList<>();
        Set<String> uniqueSampleNames = new TreeSet<>();
        ArrayList<int[]> callableBases = new ArrayList<>();
        Set<String> unannotatedGenes = new TreeSet<>();
        reader = new BufferedReader(new FileReader(sampleWigFile));
        line = reader.readLine();
        int lineNumber = 1;
        try {
            while(line != null) {
                if(line.trim().isEmpty() == false && line.startsWith("#") == false) {
                    String[] tokens = line.split("\t");
                    if(tokens.length != 2) {
                        throw new RuntimeException("Error occurred while reading sample_wig_files: " + sampleWigFile +
                                "\nat Line " + lineNumber + ": " + line + "\nError details: line must be a two-column tab-delimited string" + 
                                "\nColumn 1: Sample name or id       Column 2: wig file path");
                    }
                    if(uniqueSampleNames.add(tokens[0]) == true) {
                        sampleNames.add(tokens[0]);
                        System.out.println("Processing sample wig file: " + tokens[1] + " at line " + lineNumber + " in sample_wig_files");
                        Map<String, Integer> sampleCallableInfo = processWigSampleFile(tokens[1], geneInfo);
                        // save results of processing this sample
                        int[] result = new int[geneList.size()+1];
                        result[0] = sampleCallableInfo.get(WHOLE_EXOME_GENE_NAME);
                        for(int i=0; i<geneList.size(); i++) {
                            String gene = geneList.get(i);
                            try {
                                result[i+1] = sampleCallableInfo.get(gene);
                            } catch(NullPointerException npe) {
                                result[i+1] = -1;
                                if(unannotatedGenes.add(gene) == true) {
                                    System.out.println("Warning: " + gene + " has no annotation");
                                }
                            }
                        }
                        callableBases.add(result);
                    } else {
                        // duplicated sample name (id)
                        throw new RuntimeException("Error occurred while reading sample_wig_files: " + sampleWigFile +
                                "\nat Line " + lineNumber + ": " + line + "\nError details: duplicated sample name: " + tokens[0]);
                    }
                }
                lineNumber++;
                line = reader.readLine();
            }
        } catch(Exception ex) {
            reader.close();
            throw ex;
        }
        reader.close();
        System.out.println("Finished processing sample wig files");
        System.out.println();
        
        // Step 4: writing results
        System.out.println("Step 4: writing results to output file: " + outputFile);
        FileWriter writer = new FileWriter(outputFile);
        // Writing header line
        writer.write("Sample name\t" + WHOLE_EXOME_GENE_NAME);
        for(String gene: geneList) {
            writer.write("\t" + gene);
        }
        writer.write(LINE_SEPARATOR);
        // Writing sample by sample
        for(int i=0; i<sampleNames.size(); i++) {
            writer.write(sampleNames.get(i));
            int[] sampleResults = callableBases.get(i);
            for(int j=0; j<sampleResults.length; j++) {
                writer.write("\t" + sampleResults[j]);
            }
            writer.write(LINE_SEPARATOR);
        }
        writer.close();
        System.out.println("Finished writing results to output file: " + outputFile);
    }
    
    // Step 5 : sum up the final result
    
    
    
    public static Map<String, GeneLocationInfo> getGenesInfo(String gtfFile, ArrayList<String> geneList) throws IOException {
        Map<String, GeneLocationInfo> result = new HashMap<>();
        Set<String> genes = null;
        if(geneList == null || geneList.isEmpty() == true) {
            // genes = null; do nothing, genes is already null
        } else {
            genes = new TreeSet<>();
            for(String gene:geneList) {
                genes.add(gene);
            }
        }
        BufferedReader reader = new BufferedReader(new FileReader(gtfFile));
        String line = reader.readLine();
        while(line != null) {
            if(line.startsWith("#") == false) {
                String[] tokens = line.split("\t");
                if(tokens[2].equals("gene") == true) {
                    String info = tokens[8];
                    int nameIndex = info.indexOf("gene_id");
                    if(nameIndex > -1) {
                        // gene name is known
                        int nameStart = info.indexOf("\"", nameIndex);
                        int nameEnd = info.indexOf("\"", nameStart+1);
                        String geneName = info.substring(nameStart+1, nameEnd);
                        if(genes == null || genes.contains(geneName) == true) {
                            int geneStart = Integer.parseInt(tokens[3]);
                            int geneEnd = Integer.parseInt(tokens[4]);
                            result.put(geneName, new GeneLocationInfo(tokens[0], geneStart, geneEnd));
                        }
                    }
                }
            }
            line = reader.readLine();
        }
        reader.close();
        return result;
    }
    
    public static Map<String, Integer> processWigSampleFile(String wigFile, Map<String, GeneLocationInfo> geneInfo) throws IOException {
        if(geneInfo == null || geneInfo.isEmpty() == true) {
            throw new IllegalArgumentException("Empty gene list passed to processWigSampleFile");
        }
        Set<String> unstartedGenes = new TreeSet<>();
        for(String geneName: geneInfo.keySet()) {
            unstartedGenes.add(geneName);
        }
        Set<String> startedGenesInsideCurrentChr = new TreeSet<>();
        Set<String> unstartedGenesInsideCurrentChr = new TreeSet<>();
        int wholeGenomeCallable = 0;
        String currentChr = "";
        int currentPos = 0;
        int currentStep = 0;
        Map<String, Integer> result = new TreeMap<>();
        
        BufferedReader reader = new BufferedReader(new FileReader(wigFile));
        String line = reader.readLine();
        Set<String> tempGenes = new HashSet<>();
        while(line != null) {
            if(line.equals("1") == true || line.equals("0") == true) {
                // first update unstarted genes if necessary
                for(String gene:unstartedGenesInsideCurrentChr) {
                    GeneLocationInfo location = geneInfo.get(gene);
                    if(location.getStart() > currentPos) {
                        // still not started, do nothing
                    } else if(location.getEnd() >= currentPos) {
                        // this gene just started, add it to started and temp genes
                        startedGenesInsideCurrentChr.add(gene);
                        //tempGenes.add(gene);
                    } else {
                        // this gene has to be skipped (its position hasn't appeared or skipped in the wig file)
                        // add it to temp genes so it can be removed from started genes later
                        tempGenes.add(gene);
                    }
                }
                unstartedGenesInsideCurrentChr.removeAll(tempGenes);
                // empty tempGenes so it can be used for updating started genes below
                tempGenes.clear();
                // now update started genes if necessary
                for(String gene: startedGenesInsideCurrentChr) {
                    GeneLocationInfo location = geneInfo.get(gene);
                    if(location.getEnd() < currentPos) {
                        tempGenes.add(gene);
                    }
                }
                startedGenesInsideCurrentChr.removeAll(tempGenes);
                // empty tempGenes so it can be used in next iteration
                tempGenes.clear();
                // finally, update callable bases if necessary
                if(line.equals("1") == true) {
                    wholeGenomeCallable++;
                    for(String gene:startedGenesInsideCurrentChr) {
                        if(result.containsKey(gene) == true) {
                            result.put(gene, (result.get(gene) + 1));
                        } else {
                            result.put(gene, 1);
                        }
                    }
                } // else this base is not callable, no need to update callable bases
                currentPos += currentStep;
            } else if(line.contains("step") == true) {
                // this is a skip to a new position line
                String[] tokens = line.split(" ");
                String chr = tokens[1].substring(tokens[1].indexOf("=")+1);
                currentPos = Integer.parseInt(tokens[2].substring(tokens[2].indexOf("=")+1));
                currentStep = Integer.parseInt(tokens[3].substring(tokens[3].indexOf("=")+1));
                if(chr.equals(currentChr) == false) {
                    startedGenesInsideCurrentChr.clear();
                    unstartedGenesInsideCurrentChr.clear();
                    for(String gene:unstartedGenes) {
                        GeneLocationInfo location = geneInfo.get(gene);
                        if(location.getChromosome().equals(chr) == true) {
                            unstartedGenesInsideCurrentChr.add(gene);
                        }
                    }
                    // remove current chromosome's genes from unstarted genes
                    unstartedGenes.removeAll(unstartedGenesInsideCurrentChr);
                    currentChr = chr;
                }
            }
            line = reader.readLine();
        }
        reader.close();
        // At the end, all genes not in the result map have zero callable bases
        for(String geneName: geneInfo.keySet()) {
            if(result.containsKey(geneName) == false) {
                result.put(geneName, 0);
            }
        }
        result.put(WHOLE_EXOME_GENE_NAME, wholeGenomeCallable);
        return result;
    }
    
    public static void printUsage() {
        System.out.println("#######################################################################################################################################");
        System.out.println("## Usage: java -cp Pancancer.jar pancancer.GetCallableCounts annotation_file gene_list_file sample_wig_files output_file             ##");
        System.out.println("## Commad line arguments:                                                                                                            ##");
        System.out.println("##    annotation_file:  a gtf format file to extract gene locations for the gene list                                                ##");
        System.out.println("##    gene_list_file:   a file with gene names (must be names not ids) for which callable based to be counted                        ##");
        System.out.println("##                      Format: one gene name per line, comment lines, if any, start with # symbol                                   ##");
        System.out.println("##    sample_wig_files: a file with a list of samples and their wig file paths output by Mutect containing callable base information ##");
        System.out.println("##                      Format: a two-column tab-delimited file, comment lines, if any, start with # symbol                          ##");
        System.out.println("##                      Column 1: Sample name or id       Column 2: wig file path                                                    ##");
        System.out.println("##    output_file: a two-column tab-delimited file of the program results as follows:                                                ##");
        System.out.println("##                      Column 1: gene name (or whole exome)    Column 2: callable base count                                        ##");
        System.out.println("#######################################################################################################################################");
    }
    
    public static class GeneLocationInfo {
        private String chromosome;
        private int start;
        private int end;
        
        public GeneLocationInfo(String chromosome, int start, int end) {
            this.chromosome = chromosome;
            this.start = start;
            this.end = end;
        }

        public String getChromosome() {
            return chromosome;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }
    
    public static final String WHOLE_EXOME_GENE_NAME = "Whole Exome";
    public static final String LINE_SEPARATOR = System.lineSeparator();
}
