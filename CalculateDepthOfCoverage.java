import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Collections;

// The java version of the target_reads.py scripts.
// The script include BWA analysis result
// Each sample takes around 15 min vs 2 hr in python script
public class CalculateDepthOfCoverage {

	public static void main(String[] args) throws IOException {

		String ExonInterval = "/work/szlab/Lab_shared_PanCancer/source/Uniq_Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
		//"C:\\Users\\abc73_000\\Desktop\\CoverageTest\\Uniq_Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
		// "C:\\Users\\abc73_000\\Desktop\\CoverageTest\\Uniq_Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
		// "/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
		String DepthOfCoverage = args[0];
		//"C:\\Users\\abc73_000\\Desktop\\CoverageTest\\chr1-10.txt";
		// "C:\\Users\\abc73_000\\Desktop\\CoverageTest\\testCoverage.txt";
		// args[0];
		File file_output = new File(args[1]);
		//"C:\\Users\\abc73_000\\Desktop\\CoverageTest\\testFile.txt"
		// Create Total CDS-Interval-Dictionary
		Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict = createCDSIntervalDict(ExonInterval);

		Set<String> Total_chrom_list = Total_interval_dict.keySet();

		for (String ele : Total_chrom_list) {
			Collections.sort(Total_interval_dict.get(ele));
			// System.out.println(Total_interval_dict.get(ele).size());

		}

		calculateMeanDepth(DepthOfCoverage, file_output, Total_interval_dict);

	}

	/**
	 * @param ExonInterval
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws NumberFormatException
	 */
	private static Map<String, ArrayList<ExonLocationInfo>> createCDSIntervalDict(String ExonInterval)
			throws FileNotFoundException, IOException, NumberFormatException {
		Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict = new TreeMap<>();
		BufferedReader objReader = null;
		objReader = new BufferedReader(new FileReader(ExonInterval));
		String strCurrentLine;
		while ((strCurrentLine = objReader.readLine()) != null) {
			ArrayList<ExonLocationInfo> ExonLocation = new ArrayList<ExonLocationInfo>();

			String[] tokens = strCurrentLine.split(":");
			String chrom = tokens[0];
			int start = Integer.parseInt(tokens[1].split("-")[0]);
			int end = Integer.parseInt(tokens[1].split("-")[1]);

			ExonLocation.add(new ExonLocationInfo(start, end));
			if (Total_interval_dict.containsKey(chrom) == false) {

				Total_interval_dict.put(chrom, ExonLocation);
			} else {
				Total_interval_dict.get(chrom).add(new ExonLocationInfo(start, end));

			}

		}
		objReader.close();
		return Total_interval_dict;
	}

	public static void calculateMeanDepth(String DepthFile, File file_output,
			Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict) throws NumberFormatException, IOException {
		ArrayList<Integer> ChunkCoverage = new ArrayList<Integer>();
		BufferedWriter bw = null;
		FileWriter fw = new FileWriter(file_output);
		bw = new BufferedWriter(fw);
		int oldExonIndex = 0;
		String oldChrom = "chr1";
		try (FileReader fileReader = new FileReader(DepthFile)) {
			try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {

				String strCurrentLine;

				while ((strCurrentLine = bufferedReader.readLine()) != null) {
					if (strCurrentLine.startsWith("L")) {
						;
					} else {
						ArrayList<String> summary = analyzeLine(strCurrentLine, Total_interval_dict);
						if (!summary.isEmpty()) {
							String chrom = summary.get(0);
							int coverage = Integer.valueOf(summary.get(1));
							int newExonIndex = Integer.valueOf(summary.get(2));
							
							if (oldChrom.equals(chrom)== true && (oldExonIndex == newExonIndex)) {
								ChunkCoverage.add(coverage);
								
								
							} else if (oldChrom.equals(chrom)== false && !(oldExonIndex == newExonIndex) ) {

								double meanCoverage = calCulatedMean(ChunkCoverage);
								int startpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).start;
								int endpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).end;
								bw.write(oldChrom + ":" + String.valueOf(startpos) + "-" + String.valueOf(endpos) + "\t"
										+ String.valueOf(meanCoverage) + "\n");
								oldChrom = chrom;
								oldExonIndex = newExonIndex;
								ChunkCoverage.clear();
								ChunkCoverage.add(coverage);

							} else if (oldChrom.equals(chrom)== true && !(oldExonIndex == newExonIndex)) {
								double meanCoverage = calCulatedMean(ChunkCoverage);
								int startpos = Total_interval_dict.get(chrom).get(oldExonIndex).start;
								int endpos = Total_interval_dict.get(chrom).get(oldExonIndex).end;
								oldExonIndex = newExonIndex;
								bw.write(chrom + ":" + String.valueOf(startpos) + "-" + String.valueOf(endpos) + "\t"
										+ String.valueOf(meanCoverage) + "\n");
								ChunkCoverage.clear();
								ChunkCoverage.add(coverage);

							}
							else if (oldChrom.equals(chrom)== false && (oldExonIndex == newExonIndex)) {
								double meanCoverage = calCulatedMean(ChunkCoverage);
								int startpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).start;
								int endpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).end;
								bw.write(oldChrom + ":" + String.valueOf(startpos) + "-" + String.valueOf(endpos) + "\t"
										+ String.valueOf(meanCoverage) + "\n");
								oldChrom = chrom;
								oldExonIndex = newExonIndex;
								ChunkCoverage.clear();
								ChunkCoverage.add(coverage);
							}

						}else if (summary.isEmpty() && ChunkCoverage.size()> 0){
							double meanCoverage = calCulatedMean(ChunkCoverage);
							int startpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).start;
							int endpos = Total_interval_dict.get(oldChrom).get(oldExonIndex).end;
							bw.write(oldChrom + ":" + String.valueOf(startpos) + "-" + String.valueOf(endpos) + "\t"
									+ String.valueOf(meanCoverage) + "\n");
							ChunkCoverage.clear();
						}

					}

				}
				bw.close();
				bufferedReader.close();
			}

		}
	}

	public static double calCulatedMean(ArrayList<Integer> chunkArrayCoverage) {
		double sum = 0;
		for (int i = 0; i < chunkArrayCoverage.size(); i++) {
			sum += chunkArrayCoverage.get(i);
		}

		return (sum / chunkArrayCoverage.size());

	}

	// Analyzed for each line
	public static ArrayList<String> analyzeLine(String strCurrentLine,
			Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict) {

		ArrayList<String> CurrentLineSummary = new ArrayList<String>();
		String[] eachLine = strCurrentLine.split("\t");
		String chrom = eachLine[0].split(":")[0];
		int pos = Integer.valueOf(eachLine[0].split(":")[1]);
		int coverage = Integer.valueOf(eachLine[1]);

		ArrayList<ExonLocationInfo> Total_exome_loc = Total_interval_dict.get(chrom);

		int location_status = binarySearch(pos, Total_exome_loc, 0, (Total_exome_loc.size() - 1));

		if (location_status != -2) {
			// transcript_list.add(reads_name);
			CurrentLineSummary.add(chrom);
			CurrentLineSummary.add(String.valueOf(coverage));
			CurrentLineSummary.add(String.valueOf(location_status));
		}
		return CurrentLineSummary;
	}

	public static int withinRegion(int reads_position, int exom_start, int exom_end) {

		int value = 0;

		// means the read is greater
		if (reads_position > exom_end) {
			value = 0;
		}
		// means the read is smaller
		else if (reads_position < exom_start) {
			value = -1;
		}
		// means the read is in the region
		else {
			value = 1;
		}
		return value;

	}

	public static int binarySearch(int reads_position, ArrayList<ExonLocationInfo> arr, int left, int right) {

		// Check base case
		if (right >= left) {

			int mid = (left + right) / 2;
			int exon_start_value = arr.get(mid).getstart();
			int exon_end_value = arr.get(mid).getend();
			// If element is present at the middle itself
			int region_status = withinRegion(reads_position, exon_start_value, exon_end_value);
			if (region_status == 1) {
				return mid;
			}
			// If element is smaller than mid, then it can only be present in left subarray
			else if (region_status == -1) {
				return binarySearch(reads_position, arr, left, mid - 1);
			}

			// Else the element can only be present in right subarray
			else if (region_status == 0) {
				return binarySearch(reads_position, arr, mid + 1, right);
			}
		}

		// Element is not present in the array
		return -2;

	}
}

class ExonLocationInfo implements Comparable<ExonLocationInfo> {

	int start;
	int end;

	public ExonLocationInfo(int start, int end) {
		this.start = start;
		this.end = end;
	}

	public int getstart() {
		return start;
	}

	public int getend() {
		return end;
	}

	public int compareTo(ExonLocationInfo exon) {
		// use start point to sort
		return this.start - exon.start;
	}

}
