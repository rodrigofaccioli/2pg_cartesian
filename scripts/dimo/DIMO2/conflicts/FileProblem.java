package conflicts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class FileProblem extends Problem {		
	/* number of objectives is os_dim */
	private int n; 		        // number of points
	private double[][] points;  // objective values of the n points
	
	public FileProblem(String filename) {
		this.points = this.getPointsFromFile(filename);		
	}
	
	/**
	 *  returns the objective values for the genotype
	 */
	public double[] getFitness(Object genotype) {
		int id = ((Individual)genotype).getID();
		for (int i=0; i<points.length; i++) {
			if (points[i][0] == id) {
				double[] p = new double[points[i].length-1];
				for (int j=0; j<points[i].length-1; j++) {
					p[j] = points[i][j+1];
				}
				return p;
			}
		}		
		return null;
	}
	
	/**
	 *  return the matrix with id's in first column and the objective values for the
	 * points in the other columns
	 */
	public double[][] getPoints() {
		return this.points;
	}
	
	public int getNumberOfDifferentPoints() {
		return this.n;
	}
	
	/* returns true iff str contains a single integer */
	public static boolean isParsableToInteger(String str) {
		try {
			Integer.parseInt(str);
			return true;
		}
		catch(NumberFormatException nfe) {
			return false;
		}
	}
	
	/* splits line along the tabs if at least one tab is contained;
	 * otherwise, splits line along blanks
	 */
	public static String[] splitLine(String line) {
		if (line.contains("\t")) {
			return line.split("\t");
		} else {
			return line.split(" ");
		}
	}
	
	
	/**
	 * Reads points from file filename where the format can be one of the following.
	 * 
	 * format standard (n: #points, k: #objectives)
	 * --------------------------------------------
	 * n= 3
	 * k= 4
	 * 1 2.4 4.2 6.3 1.23
	 * 3 2.1 2.2 5.76 2.0
	 * 7 4.9 7.2 1.2 3.5
	 * EOF
	 * 
	 *  
	 * format PISA (first line: #doubles to read=n*(k+1))
	 * --------------------------------------------------
	 * 15
	 * 1 2.4 4.2 6.3 1.23
	 * 3 2.1 2.2 5.76 2.0
	 * 7 4.9 7.2 1.2 3.5
	 * END
	 * 
	 * 
	 */
	private double[][] getPointsFromFile(String filename) {
		double[][] p = new double[2][2];
		
		boolean pisaFormat = false; // true iff PISA format is detected
		
		try {
			String line;
	
			File inputFile = new File(filename);
			FileReader inputStream = new FileReader(inputFile);
			BufferedReader input = new BufferedReader(inputStream);
			
			for (int i=0; i<2; i++) {
				line = input.readLine();
				
				String[] lineSegment = splitLine(line);
				if (line.contains("k=")) {
					this.os_dim = new Integer(lineSegment[1]).intValue();
				} else if (line.contains("n=")) {
					this.n = new Integer(lineSegment[1]).intValue();
				} else if (isParsableToInteger(lineSegment[0])) {
					/* only preliminary assignment before to read first individual: */
					this.n = Integer.parseInt(lineSegment[0]);
					pisaFormat = true;
					break;
				}
			}
			if (pisaFormat) {
				/* read first line to get number of objectives */
				line = input.readLine();
				String[] lineSegment = splitLine(line);
				this.os_dim = lineSegment.length - 1;
				this.n = this.n/(this.os_dim+1); // correct number of points now
				p = new double[this.n][this.os_dim+1];
				
				/* put first individual already into p */
				for (int j=0; j<this.os_dim+1; j++) {
					p[0][j] = Double.valueOf(lineSegment[j]).doubleValue();
				}
				
				int i=1; // number of points read
				
				while ((line = input.readLine()) != null) {
					lineSegment = splitLine(line);
					if (lineSegment.length == (this.os_dim+1)){
						for (int j=0; j<this.os_dim+1; j++) {
							p[i][j] = Double.valueOf(lineSegment[j]).doubleValue();
						}
						i++;
					}
				}
				
			} else {
				p = new double[this.n][this.os_dim+1];
				int i=0; // number of points read
			
				while ((line = input.readLine()) != null) {
					String[] lineSegment = splitLine(line);
					if (lineSegment.length == (this.os_dim+1)){
						for (int j=0; j<this.os_dim+1; j++) {
							p[i][j] = Double.valueOf(lineSegment[j]).doubleValue();
						}
						i++;
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return p;
	}

}
