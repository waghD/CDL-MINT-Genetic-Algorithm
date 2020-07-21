package main;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.sql.Timestamp;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;

//import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

import harmony.HarmonyMemory;
import harmony.HarmonyParameters;
import harmony.HarmonyResult;
import harmony.HarmonySearch;

import java.util.HashMap;
import output.Printer;
import timeSeries.TimeSeriesDatabase;

public class Main_harmony {

	// Set Output directory for Test
	private static final String OUTPUT_DIRECTORY = "../harmonyresult";
	private static final String TEST_NAME = "test_refactor";

	public static void main(String[] args) {
		
		// SetUp files for Test case
		setUpDatabase("./lib/Daten_1560.csv", false, 0);
		Evaluation eval = new Evaluation("./lib/realStates_1560.csv");

		// Define Axis used by Test case
		AxisStream[] axisArr = { AxisStream.GP, AxisStream.MAP, AxisStream.BP, AxisStream.SAP, AxisStream.WP };
		
		List<AxisStream> axisList = new ArrayList<AxisStream>(Arrays.asList(axisArr));

		eval.setUpRealDataStream(axisList);
		
		// Create output directory
		String outputDir = OUTPUT_DIRECTORY + "/" + TEST_NAME;
		new File(outputDir).mkdirs();

		PrintStream out = null;
		
		// Define parameter combinations for test
		List<Double> accList = new ArrayList<Double>(Arrays.asList(0.9));
		List<Double> adjList = new ArrayList<Double>(Arrays.asList(0.5));
		List<Integer> sizeList = new ArrayList<Integer>(Arrays.asList(50));
		List<Double> bandwidthList = new ArrayList<Double>(Arrays.asList(0.4));
		
		DecimalFormat df = new DecimalFormat("#.###");
		df.setRoundingMode(RoundingMode.CEILING);

		// Run Harmony Search Algorithm with all Parameter combinations
		for (double acc : accList) {
			for (double adj : adjList) {
				for (int size : sizeList) {
					for (double bandwidth : bandwidthList) {
						
						PrintStream metaFile = null;
						
						String fileBase = outputDir + "/" + TEST_NAME + "_" + acc + "_" + adj + "_" + size + "_" + bandwidth;
						try { 
							metaFile = new PrintStream(new FileOutputStream(fileBase + "_main" + ".txt", true), true); 
							System.setOut(metaFile); 
						} catch (IOException e) { 
							System.err.print(e.getMessage());
							e.printStackTrace();
						}

						// states NOT to use for evaluation
						List<String> statesToNotEvaluateList = new ArrayList<String>();
						List<HarmonyResult> resultList = new ArrayList<HarmonyResult>();

						// configure Harmony Search
						HarmonyParameters hpa = new HarmonyParameters(adj, bandwidth, acc, size, axisList);
						hpa.setPrintNewSolutions(false);
						hpa.setPrintMemorySwaps(false);
						hpa.setNrOfIterations(500);
						hpa.setStopOnOptimum(false);
						hpa.setStatesToNotEvaluate(statesToNotEvaluateList);
						hpa.setMinimizeResult(true);
						hpa.setLowerSearchBorder(0.0);
						hpa.setUpperSearchBorder(0.4);

						// Print Harmony Parameters to main file
						Printer.printHeader("Harmony parameters:");
						System.out.println(hpa);

						// Print Headers for output csv
						System.out.println("Iteration;AvgFMeasure;AvgPrecision;AvgRecall;IterationBestFMeasure;IterationMinimization");
						
						// number of iterations for average calculation
						for (int i = 0; i < 5; i++) {
							// Set output stream to file for this run
							try { 
								out = new PrintStream(new FileOutputStream(fileBase + "_" + i + ".txt", true), true); 
								System.setOut(out); 
							} catch (IOException e) { 
								System.err.print(e.getMessage());
								e.printStackTrace();
							}
							
							// Run Harmony Search Algorithm
							HarmonyResult result = runHarmonySearch(hpa);
							out.close();
							
							// Print averages to main file
							System.setOut(metaFile);
							System.out.print(i + ";" + result.getAvgBestFMeasure() + ";" + result.getAvgBestPrecision() + ";" + result.getAvgBestRecall() +";"+ result.getNrOfIterationsForBestFMeasure() + ";" + result.getNrOfIterationsForBestMinimizedRange() + "\n");
							resultList.add(result);
						}

						System.setOut(metaFile); 

						// Collect results for all runs
						List<Integer> iterationsFMeasureList = new ArrayList<Integer>();
						List<Integer> iterationsMinRangeList = new ArrayList<Integer>();
						List<Double> absOffsetPercDiff = new ArrayList<Double>();
						List<Double> timeTilOptList = new ArrayList<Double>();
						List<Double> avgPrecisionList = new ArrayList<Double>();
						List<Double> avgRecallList = new ArrayList<Double>();
						List<Double> avgFMeasureList = new ArrayList<Double>();
						List<Double> timeOverallList = new ArrayList<Double>();
						int repetition = 0;
						for (HarmonyResult res : resultList) {
							System.out.println("Iterations til best f-measure (Exec. " + repetition + "): " + res.getNrOfIterationsForBestFMeasure() + " ("
									+ res.getRuntimeTilOptimumFound() + "s)"
							);
							
							iterationsFMeasureList.add(res.getNrOfIterationsForBestFMeasure());
							iterationsMinRangeList.add(res.getNrOfIterationsForBestMinimizedRange());
							absOffsetPercDiff.add(res.getAbsOffsetBestMinimized()/res.getAbsOffsetBestInitial());
							timeTilOptList.add(res.getRuntimeTilOptimumFound());
							avgPrecisionList.add(res.getAvgBestPrecision());
							avgRecallList.add(res.getAvgBestRecall());
							avgFMeasureList.add(res.getAvgBestFMeasure());
							timeOverallList.add(res.getRuntimeIterations());
							repetition++;
						}

						// Calculate averages over all runs
						double avgIterationsFMeasure = iterationsFMeasureList.stream().mapToInt(x -> x).average().orElse(-1);
						double avgIterationsMinRange = iterationsMinRangeList.stream().mapToInt(x -> x).average().orElse(-1);
						double avgAbsOffsetPercDiff = absOffsetPercDiff.stream().mapToDouble(x -> x).average().orElse(-1);

						double avgTimeTilOpt = timeTilOptList.stream().mapToDouble(x -> x).average().orElse(-1);
						double avgTimeOverall = timeOverallList.stream().mapToDouble(x -> x).average().orElse(-1);
						double avgPrec = avgPrecisionList.stream().mapToDouble(x -> x).average().orElse(-1);
						double avgRec = avgRecallList.stream().mapToDouble(x -> x).average().orElse(-1);
						double avgFMeasure = avgFMeasureList.stream().mapToDouble(x -> x).average().orElse(-1);

						// Print averages to main file
						System.out.println(Printer.div);
						Printer.printHeader("After all repetitions:");
						System.out.println("Iterations til best f-measure (avg): " + avgIterationsFMeasure);
						System.out.println("Iterations til best minimized range (avg): " + avgIterationsMinRange);
						System.out.println("Abs Offset Difference (between best f-measure and offset optimized best f-measure), (avg): " + ((1.0-avgAbsOffsetPercDiff)*100.0) + "%");
						System.out.println("Time til Best found (avg): " + avgTimeTilOpt);
						System.out.println("Time overall (avg): " + avgTimeOverall);
						System.out.println("Avg Precision: " + avgPrec);
						System.out.println("Avg Recall: " + avgRec);
						System.out.println("Avg F-measure: " + avgFMeasure);

						if(out != null) {
							out.close();
						}
						if(metaFile != null) {
							metaFile.close();
						}
					}
				}
			}
		}
	}

	/**
	 * Execute harmony search with given parameters and prints some information
	 * (influenced by constants in main).
	 * 
	 * @param hpa .. Configuration Object of Harmony Search
	 * 
	 * @return HarmonyResult
	 */
	static HarmonyResult runHarmonySearch(HarmonyParameters hpa) {
		// Initialize HarmonyMemory with HarmonyParameters
		HarmonyMemory memory = new HarmonyMemory(hpa);

		// Initialize Harmony Search with Harmony Memory and Parameters
		HarmonySearch search = new HarmonySearch(memory, hpa);

		// Execute harmony search: If parameter "stopIfOptimumFound" is true,
		// max number of iterations is still respected (here: max 300 iterations)
		HarmonyResult hs = search.execHarmonySearch();

		// Print results of this Harmony Search run
		Printer.printHeader("HARMONY RESULTS");

		System.out.println(hs);
		Printer.printHeader("BEST");
		memory.print(memory.findBestEvalResult());

		return hs;
	}


	private static void setUpDatabase(String filenameData, boolean longRun, int timespan) {
		TimeSeriesDatabase db = new TimeSeriesDatabase();
		db.setUpData(filenameData, longRun, timespan);
	}
}
