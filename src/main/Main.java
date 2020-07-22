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

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static io.jenetics.engine.EvolutionResult.toBestPhenotype;
import static io.jenetics.engine.Limits.bySteadyFitness;

import io.jenetics.Chromosome;
import io.jenetics.DoubleChromosome;
import io.jenetics.DoubleGene;
import io.jenetics.Genotype;
import io.jenetics.MeanAlterer;
import io.jenetics.Mutator;
import io.jenetics.Optimize;
import io.jenetics.Phenotype;
import io.jenetics.RouletteWheelSelector;
import io.jenetics.TournamentSelector;
import io.jenetics.engine.Codecs;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.engine.EvolutionStatistics;
import io.jenetics.stat.MinMax;
import io.jenetics.util.DoubleRange;
//import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;
import io.jenetics.util.ISeq;
import harmony.HarmonyMemory;
import harmony.HarmonyParameters;
import harmony.HarmonyResult;
import harmony.HarmonySearch;

import java.util.HashMap;
import output.Printer;
import timeSeries.TimeSeriesDatabase;

public class Main {

	/*
	 * private static final String OUTPUT_DIRECTORY = "../harmonyresult"; private
	 * static final String TEST_NAME = "test_refactor";
	 */
	private static final AxisStream[] axisArr = { AxisStream.GP, AxisStream.MAP, AxisStream.BP, AxisStream.SAP,
			AxisStream.WP };
	private static List<AxisStream> axisList = null;
	private static final List<String> statesToNotEvaluateList = new ArrayList<String>();
	
	// Test Output
	private static final boolean OUTPUT_TO_FILE = true;
	private static final String OUTPUT_DIRECTORY = "../genetic_result";
	private static final String TEST_NAME = "test_output";

	// The fitness function.
	private static double fitness(final Genotype<DoubleGene> x) {

		// Allele = Jenetic's wording for the actual gene values that
		// are stored in a chromosome (in our case of type double)
		double[] alleles = x.chromosome().as(DoubleChromosome.class).toArray();

		Map<String, PropertyBoundaries> propertyMap = getSensorOffsetMap(alleles);

		List<EvaluationResult> resultList = performStateDetection(propertyMap);
		double fMeasure = 0.0;
		for (EvaluationResult res : resultList) {
			fMeasure += res.getfMeasure() / resultList.size();
			// newPrec += res.getPrecision() / newResult.size();
			// newRec += res.getRecall() / newResult.size();
		}
		return fMeasure;
	}

	public static void main(String[] args) {

		// Add sensors to use for detection here, derive of "AxisStream"-Class
		setUpDatabase("./lib/Daten_156.csv", false, 0);
		Evaluation eval = new Evaluation("./lib/realStates_156.csv");

		int nrOfSensors = axisArr.length;
		axisList = new ArrayList<AxisStream>(Arrays.asList(axisArr));

		eval.setUpRealDataStream(axisList);

		PrintStream outputFile = null;
		PrintStream resultFile = null;
		
		// Create output directory
		if(OUTPUT_TO_FILE) {
			String outputDir = OUTPUT_DIRECTORY + "/" + TEST_NAME;
			new File(outputDir).mkdirs();
			
			String fileBase = outputDir + "/" + TEST_NAME;
			try { 
				outputFile = new PrintStream(new FileOutputStream(fileBase + "_logs" + ".csv", true), true); 
				resultFile = new PrintStream(new FileOutputStream(fileBase + "_result" + ".txt", true), true);
				System.setOut(outputFile); 
			} catch (IOException e) { 
				System.err.print(e.getMessage());
				e.printStackTrace();
			}
		}

		// final Factory<Genotype<DoubleGene>> gtf = Genotype.of(DoubleChromosome.of())

		final Engine<DoubleGene, Double> engine = Engine
				// Create a new builder with the given fitness
				// function and chromosome.
				.builder(Main::fitness, DoubleChromosome.of(0.0, 1.0, nrOfSensors * 2))
				.offspringFraction(0.7)
				.survivorsSelector(new RouletteWheelSelector<>())
				.offspringSelector(new TournamentSelector<>())
				.populationSize(100)
				.optimize(Optimize.MAXIMUM)
				.alterers(new Mutator<>(0.1), new MeanAlterer<>(0.8))
				// Build an evolution engine with the
				// defined parameters.
				.build();

		// Create evolution statistics consumer.
		final EvolutionStatistics<Double, ?> statistics = EvolutionStatistics.ofNumber();
		
		StringBuilder headerStr = new StringBuilder();
		for(int i = 0; i < nrOfSensors; i++) {
			headerStr.append("Sensor ");
			headerStr.append(i);
			headerStr.append(" Lower;");
			headerStr.append("Sensor ");
			headerStr.append(i);
			headerStr.append(" Upper;");
		}
		
		headerStr.append("F-Measure");
		
		System.out.println(headerStr);

		final ISeq<EvolutionResult<DoubleGene, Double>> best = engine.stream()
				// Truncate the evolution stream after 7 "steady"
				// generations.
				// .limit(bySteadyFitness(5))
				// The evolution will stop after maximal 100
				// generations.
				.limit(1000)
				// Update the evaluation statistics after
				// each generation
				.peek(statistics)
				.peek((evolutionResult) -> {
					Phenotype<DoubleGene, Double> bestPheno = evolutionResult.bestPhenotype();
					Genotype<DoubleGene> genum = bestPheno.genotype();
					StringBuilder str = new StringBuilder();
					for(Chromosome<DoubleGene> chromosome : genum) {
						for(DoubleGene gene : chromosome) {
							str.append(gene.doubleValue());
							str.append(";");
						}
					}
					str.append(bestPheno.fitness());
					System.out.println(str);
				})
				// Collect (reduce) the evolution stream to
				// its best phenotype.
				// .collect(EvolutionResult.toBestEvolutionResult());
				// .flatMap(MinMax.toStrictlyDecreasing())
				.collect(ISeq.toISeq(1000));

		// .collect(ISeq.toISeq(10));
		
		if(OUTPUT_TO_FILE) {
			System.setOut(resultFile); 
		}
		
		System.out.println(statistics);
		ArrayList<EvolutionResult<DoubleGene, Double>> arrayList = new ArrayList<>(best.asList());
		double fittestVal = 0.0;
		for (int i = 0; i < arrayList.size(); i++) {
			EvolutionResult<DoubleGene, Double> curResult = arrayList.get(i);
			List<Phenotype<DoubleGene, Double>> solutions = curResult.population().asList();
			for (int j = 0; j < solutions.size(); j++) {
				Phenotype<DoubleGene, Double> curSol = solutions.get(j);
				if (curSol.fitness() > fittestVal) {
					System.out.println(i + ":");
					System.out.println(curSol);
					fittestVal = curSol.fitness();
				}

			}
		}
		System.out.println(fittestVal);

		/*
		 * for(EvolutionResult<DoubleGene, Double> result: best) {
		 * System.out.println(result.population());
		 * System.out.println(result.population().size()); }
		 */
		if(outputFile != null) {
			outputFile.close();
		}
		if(resultFile != null) {
			outputFile.close();
		}
	}

	private static List<EvaluationResult> performStateDetection(Map<String, PropertyBoundaries> sensorOffsetMap) {
		Evaluation eval = Evaluation.instance;
		// evaluate new solution
		List<EvaluationResult> evalResults = eval.evaluate(TestData.setUpDataStream(axisList), sensorOffsetMap, false,
				statesToNotEvaluateList);
		return evalResults;
	}

	private static Map<String, PropertyBoundaries> getSensorOffsetMap(double[] alleles) {

		Map<String, PropertyBoundaries> sensorOffsetMap = new HashMap<String, PropertyBoundaries>();

		if (axisList.contains(AxisStream.BP)) {
			sensorOffsetMap.put(AxisStream.BP.getAxisName(), new PropertyBoundaries(alleles[0], alleles[1]));
		}

		if (axisList.contains(AxisStream.GP)) {
			sensorOffsetMap.put(AxisStream.GP.getAxisName(), new PropertyBoundaries(alleles[2], alleles[3]));

		}

		if (axisList.contains(AxisStream.MAP)) {
			sensorOffsetMap.put(AxisStream.MAP.getAxisName(), new PropertyBoundaries(alleles[4], alleles[5]));

		}

		if (axisList.contains(AxisStream.SAP)) {
			sensorOffsetMap.put(AxisStream.SAP.getAxisName(), new PropertyBoundaries(alleles[6], alleles[7]));

		}

		if (axisList.contains(AxisStream.WP)) {
			sensorOffsetMap.put(AxisStream.WP.getAxisName(), new PropertyBoundaries(alleles[8], alleles[9]));

		}
		return sensorOffsetMap;
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
