
package conflicts;
/*
 * Takes files coming from a PISA module (e.g. dtlz_output.log) and transfers
 * them to the format of input data files for FirstExperiments, i.e. it
 *    - adds the two lines with 'k=' and 'n=' at the beginning of the file
 *    - adds an EOF at the end of the file
 * 
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

public class FilePreparation {
	private int k = 1;
	private int n = 1;

	private void start(String filename, String outputfilename) {
		getKandN(filename);
		transform(filename, outputfilename);
	}

	private void getKandN(String filename) {
		try {
			String line;
			int numberoflines = 0;
			

			File inputFile = new File(filename);
			FileReader inputStream = new FileReader(inputFile);
			BufferedReader input = new BufferedReader(inputStream);
			String[] lineSegment = null;

			while ((line = input.readLine()) != null) {
				lineSegment = line.split(" ");
				numberoflines++;
			}
			if (lineSegment != null) {
				this.k = lineSegment.length - 2;
			}
			this.n = numberoflines;
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void transform(String filename, String outputfilename) {
		Vector<String> toPrint = new Vector<String>(); // to store output
		
		try {
			String line;

			File inputFile = new File(filename);
			FileReader inputStream = new FileReader(inputFile);
			BufferedReader input = new BufferedReader(inputStream);

			toPrint.add("k= " + k);
			toPrint.add("n= " + n);

			while ((line = input.readLine()) != null) {
				String[] lineSegment = line.split(" ");
				String nextLine = "";
				if (lineSegment.length > 2) {
					nextLine = nextLine + lineSegment[0] + " ";
					for (int i=1; i<lineSegment.length-1; i++) {
						nextLine = nextLine + lineSegment[i] + " ";
					}
					toPrint.add(nextLine);
				}				
			}
			toPrint.add("EOF");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* finally output everything */
		Output.print(toPrint, outputfilename);
	}
	
	/**
	 * @param args args[0]: filename
	 * 			   
	 * 	           an additional optional argument '-o outputfilename' will indicate
	 *             that all output is written to the file 'outputfilename'
	 */
	public static void main(String args[]) {
		String outputfilename = ""; // standard: output written to stdout
		if (args == null || !(args.length == 1 || args.length == 3)) {
			System.out.println("Wrong usage.");
			System.out.println();
			System.out.println("Usage:");
			System.out.println("   FilePreparation filename [-o outputfilename]");
			System.out.println();
			System.out.println("where");
			System.out.println("   filename is a file written by the PISA knapsack module");
			System.out.println();
			System.out.println("Adding '-o outputfilename' as last argument will result");
			System.out.println("   in writing the transformed knapsack output file to the");
			System.out.println("   file 'outputfilename' instead of writing to standard output.");
		} else {
			FilePreparation fp = new FilePreparation();
			if (args.length == 3) {
				outputfilename = args[2];
			}
			fp.start(args[0], outputfilename);
		}
	}

}
