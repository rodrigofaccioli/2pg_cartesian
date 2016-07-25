package conflicts;

import java.io.*;
import java.util.*;


public class Output {
	
	/*
	 * prints 'lines' to file called 'filename'
	 * 
	 * if filename is an empty string then 'lines' is printed to standard out
	 */
	public static void print(Vector<String> lines, String filename){
		if (filename.isEmpty()) {
			for (String line : lines) {
				System.out.println(line);
			}
		} else {
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(filename));
				for (String line : lines) {
					out.write(line);
					out.newLine();
				}
				out.close();
			}
			catch (IOException e) {
				System.out.println("Problems while writing output to " + filename);
			}
		}
	}

	public static void print2PG(Vector<String> lines, String filename){
                /*This method is used to integrate 2PG framework*/
                if (filename.isEmpty()) {
                        new IOException("You have to informe a file\n");
                }
                Vector<String> toPrint = new Vector<String>();
                Vector<Integer> used_objective = new Vector<Integer>();
                Vector<String> vbestElem = new Vector<String>();
                Vector<Integer> vchooseObjective = new Vector<Integer>();
                Vector<Float> verrorObjective = new Vector<Float>();

                boolean find_ = false;
                float best_error, aux_error = 0;
                int best_objective = 0, aux_objetive = 0, out_objetive = -1;
                int howMany = 0; 
                String[] aux_line, aux_line_2, aux_l;
                String aux;
                String elements = "", best_elements = "";
                best_error = -1;
                for (String line : lines) {                        
                        if ( line.startsWith("-")){
                                if (find_ == false){
                                        find_ = true;
                                }else{                                        
                                        //Saving in Vector that will be evaluated
                                        howMany = howMany + 1;
                                        vbestElem.add(best_elements);
                                        vchooseObjective.add(best_objective);
                                        verrorObjective.add(best_error);
                                        best_error = -1;
                                        best_elements = "";
                                }
                        }else{
                                
                                /* Here we will give the information */
                               aux_line = line.split("}");
                               elements = aux_line[0].replace('{',' ').trim();
                               //Split the second part                               
                               aux = aux_line[1].trim();
                               aux_line_2 = aux.split(" ",2);                                
                               aux_error = Float.parseFloat(aux_line_2[0].trim());
                               aux_objetive = Integer.parseInt(aux_line_2[1].trim());
                               if ( aux_error > best_error){
                                         best_error = aux_error;
                                         best_objective = aux_objetive;
                                         best_elements = elements;
                               }
                               
                        }
                }
                //Adding the last line
                howMany = howMany + 1;
                vbestElem.add(best_elements);
                vchooseObjective.add(best_objective);
                verrorObjective.add(best_error);
                //Evaluating vectors
                for (int j = 0; j < howMany; j++){
                         aux_l = vbestElem.get(j).split("\\s+");
                         for (int i = 0; i < aux_l.length; i++){                                                
                                 if ( (used_objective.contains(Integer.parseInt(aux_l[i])) == false) &&
                                      (Integer.parseInt(aux_l[i]) != vchooseObjective.get(j)) ){
                                              out_objetive = Integer.parseInt(aux_l[i]);
                                              break;
                                 }
                         }
                         used_objective.add(out_objetive);
                         toPrint.add(String.valueOf(verrorObjective.get(j)) + " "+ String.valueOf(out_objetive));
                }

                //Saving
                try {
			BufferedWriter out = new BufferedWriter(new FileWriter(filename));
			for (String line : toPrint) {
                                if (line.trim() != ""){
               				out.write(line);
       			        	out.newLine();
                                }
			}
			out.close();
		} catch (IOException e) {
			System.out.println("Problems while writing output to " + filename);
		}
                
	}

}
