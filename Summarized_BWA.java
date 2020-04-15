import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
 
public class Summarized_BWA {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String SamFile = args[0];
		String ID = args[1];
		String file_output = args[2];
		getSummaryBwa(SamFile,ID,file_output);
	}

	public static void getSummaryBwa(String SamFile, String ID, String file_output ) throws IOException{
		double total = 0 ;
		double unique = 0 ;//3
		double duplicate = 0 ;//3
		double Onemapped = 0 ;//5,9
		double incorrect = 0 ;//1
		double unmapped = 0 ;//13
		
		try(FileReader fileReader = new FileReader(SamFile)){
		try(BufferedReader bufferedReader= new BufferedReader(fileReader)){
		String strCurrentLine_sam;
		while ((strCurrentLine_sam = bufferedReader.readLine()) != null) {
			
		if(strCurrentLine_sam.startsWith("@")) {
			 ;
		}
		else {
			String [] file_lst =  strCurrentLine_sam.split("\t");
			total += 1 ;
            int status = Integer.parseInt(file_lst[1])%16 ;
            
            if(status == 5 || status ==9) {
            	Onemapped+=1;
            }
            else if(status ==1){
            	incorrect += 1 ;
            }
            else if(status ==13){
            	unmapped += 1;
            }
            else if (status == 3) {
		    if (file_lst[11].contains("XT:")) { 
            String status2 = file_lst[11].split(":")[2];
        if (status2.equals("U") || status2.equals("M")) {
            unique += 1 ;}
        else if (status2.equals("R")) {
        	duplicate +=1 ;}
        }
            }        
		    }
	}      
		bufferedReader.close();
		double pairs = total / 2 ;

		String unique_rate = Double.toString(unique/total) ;
		String dup_rate = Double.toString(duplicate/total) ;
		String Onemap_rate = Double.toString(Onemapped /total) ;
		String incorrect_rate = Double.toString(incorrect/total) ;
		String unmapped_rate = Double.toString(unmapped /total) ;
		BufferedWriter bw = null;
		FileWriter fw = new FileWriter(file_output);
		bw = new BufferedWriter(fw);
		//bw.write('ID\tTotal_pairs\tUniquely_mapped_rate\tRepeatedly_mapped_rate\t1read_mapped_rate\tIncorrectly_mapped_rate\tUnmapped_rate\tUniquely_mapped\tRepeatedly_mapped\t1read_mapped\tIncorrectly_mapped\tUnmapped\tTotal_reads\n')		//bw.write("File_name\t"+ "Total_reads\t"+"Total_uniq\t"+"uniq_mapped_rate\t"+"Total_read_pairs\t"+"uniq_Exonic_region\t"+"uniq_Exonic_region_paris_rates\t"+'\n');
		bw.write(ID+"\t"+pairs+"\t"+ unique_rate+"\t"+dup_rate+"\t"+Onemap_rate+"\t"+incorrect_rate+"\t"+unmapped_rate+"\t"
		+Double.toString(unique)+"\t"+Double.toString(duplicate)+"\t"+Double.toString(Onemapped)+"\t"+Double.toString(incorrect)+"\t"
		+Double.toString(unmapped));
		
		bw.close();	
	}
		
		} 
}
}