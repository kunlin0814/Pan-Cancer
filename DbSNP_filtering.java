import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class DbSNP_filtering {

	public static void main(String[] args) throws NumberFormatException, FileNotFoundException, IOException {
		long time_start = System.currentTimeMillis();
		String vcfDatabase = args[0];

		// "C:\\Users\\abc73_000\\Desktop\\DbSNP_test.vcf" ;
		String Mutect2VCF = args[1];
		// "C:\\Users\\abc73_000\\Desktop\\Consensus.sSNV.vcf";
		String file_output = args[2];
		// "C:\\Users\\abc73_000\\Desktop\\test_java.txt"

		// Create Dictionary
		Map<String, ArrayList<String>> vcfDict = createVCFDict(vcfDatabase);
		long time_Dict = System.currentTimeMillis();
		System.out.println("Create dict takes " + (time_Dict - time_start));

		BufferedWriter bw = null;
		FileWriter fw = new FileWriter(file_output);
		bw = new BufferedWriter(fw);
		// Read the VCF file
		long Compare = System.currentTimeMillis();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(Mutect2VCF), Charset.forName("UTF-8"))) {
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.startsWith("#") == false) {
					String[] tokens = line.split("\t");

					String chrom = tokens[0];
					String pos = tokens[1];
					String ref = tokens[3];
					String alt = tokens[4];
					// String status = tokens[6];
					// chrom_pos CurrentLinechrom_pos = new chrom_pos(chrom, pos);
					String CurrentLinechrom_pos = chrom + ":" + pos;
					// ref_alt CurrentLineref_alt = new ref_alt(ref, alt);
					String CurrentLineref_alt = ref + ":" + alt;
					if (vcfDict.containsKey(CurrentLinechrom_pos) == false) {
						bw.write(line + "\n");
					} else {
						ArrayList<String> Ref_alt_List = vcfDict.get(CurrentLinechrom_pos);
						if ((Ref_alt_List.contains(CurrentLineref_alt)) == false) {
							bw.write(line + "\n");
						}
					}

				}
			}
			in.close();
		} catch (IOException e) {
			// TODO: handle exception
		}
		bw.close();
		long program_end = System.currentTimeMillis();
		System.out.println("Comparing takes " + (program_end - Compare));
	}

	private static Map<String, ArrayList<String>> createVCFDict(String vcfDatabase)
			throws FileNotFoundException, IOException, NumberFormatException {
		Map<String, ArrayList<String>> VCFDict = new HashMap<>();

		// Map<String, Map<String,ArrayList<ExonLocationInfo>>> Second_dict = new

		BufferedReader objReader = null;
		objReader = new BufferedReader(new FileReader(vcfDatabase));
		String strCurrentLine;
		while ((strCurrentLine = objReader.readLine()) != null) {
			if (strCurrentLine.startsWith("#") == false) {
				String[] tokens = strCurrentLine.split("\t");
				String chrom = tokens[0];
				String pos = tokens[1];
				String ref = tokens[3];
				String alt = tokens[4];

				// chrom_pos Chrom_Pos = new chrom_pos(chrom, pos);
				// ref_alt Ref_ALT = new ref_alt(ref, alt);
				String Ref_ALT = ref + ":" + alt;
				String Chrom_Pos = chrom + ":" + pos;
				if (VCFDict.containsKey(Chrom_Pos) == false) {
					ArrayList<String> ref_alt = new ArrayList<String>();
					ref_alt.add(Ref_ALT);
					VCFDict.put(Chrom_Pos, ref_alt);
				} else {
					VCFDict.get(Chrom_Pos).add(Ref_ALT);
				}
			}
		}
		objReader.close();

		return VCFDict;
	}
}
