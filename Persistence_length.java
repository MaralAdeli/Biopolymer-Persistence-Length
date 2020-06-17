package visualize;

import inputs.Shape_inputs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.File;



import org.opensourcephysics.frames.PlotFrame;



import vector.Force;
import misc.HaosuMisc;
import misc.WriteFile;

import java.util.Collections;

public class Persistence_length {


	static int beadNumber = 0;
	static double l0 = 0.1;
	static  int length;
	static int Dim = 3;
	static int dsattempts=21;
	static  ArrayList<ArrayList<Double>> totalF = new ArrayList<ArrayList<Double>>();
	static  ArrayList<ArrayList<Double>> totalFSTDV = new ArrayList<ArrayList<Double>>();
	static ArrayList<String> time = new ArrayList<String>();


	static String[] screenCaptureTime = {};

	public static void main(String[] args) throws NumberFormatException, IOException {



		String filepath = HaosuMisc.getLoadFileName(null, "Get Trajectory File");
		String dirpath = filepath.substring(0, filepath.lastIndexOf("\\")+1);
		String filename = filepath.substring(filepath.lastIndexOf("\\")+1, filepath.length());
		Shape_inputs shape_inputs = new Shape_inputs(dirpath+"shape_inputs.txt", filepath.substring(filepath.lastIndexOf("\\")+1, filepath.length()).replaceAll("%3d", "="));
		ArrayList<ArrayList<BeadInfo>> beads = new ArrayList<ArrayList<BeadInfo>>();


		BufferedWriter[] writers = new BufferedWriter[1];
		writers[0] =  WriteFile.getWriter(filepath+"_bundleNumber.txt");

		//Generate Reader
		BufferedReader br = HaosuMisc.getReader(filepath);	//reader

		//Read Header file
		String strLine;		//reading the line...

		for (int ii = 0 ; ii < dsattempts ; ii++) {
			totalF.add(new ArrayList<Double>());
			totalFSTDV.add(new ArrayList<Double>());
		}

		boolean headerTag = true;

		while((strLine=br.readLine()) != null){
			if(strLine.contains("T=599.9")){
				while((strLine=br.readLine()) != null){
					if(strLine.contains("*****") && headerTag){	//beginning of header
						headerTag = false;
					}
					else if(strLine.contains("*****") && !headerTag){	//end of header
					}
					else if(strLine.contains("T=")){

						if(screenCaptureTime.length!=0){

							for(String s : screenCaptureTime)

								if(strLine.contains(s)){


									beads = getShape_T(strLine, br, shape_inputs, filepath, s);
								}
						}
						else{

							beads = getShape_T(strLine, br, shape_inputs, filepath, strLine);

						}



					}
				}
			}
		}

	
		double mean_correlation[]=new double[dsattempts];
		double mean_correlation_STDV[]=new double[dsattempts];
		for (int ii = 0 ; ii < dsattempts ; ii++) {
			for (int jj = 0 ; jj < totalF.get(ii).size() ; jj++) {

				mean_correlation[ii] += totalF.get(ii).get(jj);
				mean_correlation_STDV[ii] += totalFSTDV.get(ii).get(jj);

			}

		}
		double taverage[] = new double[dsattempts];
		double tstd[] = new double[dsattempts];
	
		for (int ds = 0 ; ds < dsattempts ; ds++) {
			int Ntot = totalF.get(ds).size();
		taverage[ds]=mean_correlation[ds]/(Ntot);		
				}
		
		for (int ds = 0 ; ds < dsattempts ; ds++) {
			for (int jj = 0 ; jj < totalF.get(ds).size() ; jj++) {
				int Ntot = totalF.get(ds).size();
				tstd[ds]+=((totalF.get(ds).get(jj)-taverage[ds])*(totalF.get(ds).get(jj)-taverage[ds]));

			}

		}
		
		for (int ds = 0 ; ds < dsattempts ; ds++) {
			int Ntot = totalF.get(ds).size();
			tstd[ds]=Math.sqrt(tstd[ds]/(Ntot-1)); 
			System.out.println(tstd[ds]);
		}
		PlotFrame frame = new PlotFrame("s, contour length(um)", "<t(0).t(s)>", "tangent corellation function");
		for (int i=0; i<dsattempts; i++) {

			
			frame.append(0,i*l0, taverage[i] ,0.,tstd[i]);


		}
		frame.setVisible(true);
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

	}

	private static ArrayList<ArrayList<BeadInfo>> getShape_T(String strLine, BufferedReader br, Shape_inputs shape_inputs, String filepath, String s) throws NumberFormatException, IOException{
		ArrayList<ArrayList<BeadInfo>> sBeads = new ArrayList<ArrayList<BeadInfo>>();

		beadNumber = 0;
		int NOfilament=0; 
		int count=0;


		while(!(strLine=br.readLine()).contains(";")){
			

			if(strLine.contains("RF")){


				sBeads.add(new ArrayList<BeadInfo>());

				ArrayList<Double> xcoord = new ArrayList<Double>();
				ArrayList<Double> ycoord = new ArrayList<Double>();
				ArrayList<Double> zcoord = new ArrayList<Double>();

				while(!(strLine=br.readLine()).contains("#")){

					int newbeadIndex = sBeads.size()-1;
					beadNumber++;

					double x = Double.valueOf(strLine.substring(0, strLine.indexOf(",\t")));
					double y = Double.valueOf(strLine.substring(strLine.indexOf(",\t")+2, strLine.lastIndexOf(",\t")));
					double z = Double.valueOf(strLine.substring(strLine.lastIndexOf(",\t")+2, strLine.length()));
					xcoord.add(x);
					ycoord.add(y);
					zcoord.add(z);



					sBeads.get(newbeadIndex).add(new BeadInfo(x,y,z));
				}

				length = xcoord.size();


				int Step = length-1;
				int Step_length =  length/Step;
			

				double [][] tangent = new double [Step][Dim];
				double [][] vector = new double [Step][Dim];
				double [] magnitude = new double [Step];




				for (int ii = 0; ii < Step; ii++){
					magnitude[ii] = 0;
					for (int jj=0; jj < Dim ; jj++) {
						vector [ii][jj] = 0;
						tangent [ii][jj] = 0;
					}

				}



				// defining tangent vector along each filament	


				for (int ii = 0; ii < Step; ii++){


					vector [ii][0] = (xcoord.get((ii+1)*Step_length) - xcoord.get(ii*Step_length));
					vector [ii][1] = (ycoord.get((ii+1)*Step_length) - ycoord.get(ii*Step_length));
					vector [ii][2] = (zcoord.get((ii+1)*Step_length) - zcoord.get(ii*Step_length));

					magnitude[ii] = Math.sqrt((vector[ii][0]*vector[ii][0])+(vector[ii][1]*vector[ii][1])+(vector[ii][2]*vector[ii][2]));

					tangent[ii][0] = vector[ii][0]/magnitude[ii];
					tangent[ii][1] = vector[ii][1]/magnitude[ii];
					tangent[ii][2] = vector[ii][2]/magnitude[ii]; 

				}


				int N = length-1;

				double taverage[] = new double[dsattempts];
				double tstd[] = new double[dsattempts];
				double a1[]=new double[3];
				double b1[]=new double[3];
				double c1[]=new double[3];
				double d1[]=new double[3];



				double ttotal[]=new double[dsattempts];
				double tsqtotal[]=new double[dsattempts];


				for(int ds=0;ds<dsattempts;ds++){
					for(int j=0;j<N-ds-1;j++){
						for(int i=0;i<3;i++){
							c1[i]=tangent[j+ds][i];
							d1[i]=tangent[j][i];
						}
						a1=normalizeVector(c1);
						b1=normalizeVector(d1);
						double c=vectorproduct(a1,b1);
						totalF.get(ds).add(c);

						
						totalFSTDV.get(ds).add(c*c);
			
				
				}

			
				
				}       

						}
		}

		return sBeads;


	}


	protected static double magnitude(double[] a){
		double mag = 0;
		for(int i = 0; i<3; i++)
			mag += Math.pow(a[i],2);
		return Math.sqrt(mag);
	}

	protected static double[] normalizeVector(double[] a){
		double mag = magnitude(a);  
		double[] c = new double[3];
		for(int i = 0; i<3; i++)
			c[i] = a[i]/mag;
		return c;
	}

	protected static double vectorproduct(double[] a, double[] b){
		double c=0;
		for(int i=0;i<3;i++){
			c+=a[i]*b[i];
		}
		return c;

	}

}

