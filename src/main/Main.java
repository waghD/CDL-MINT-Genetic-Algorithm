package main;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import io.jenetics.Chromosome;
import io.jenetics.DoubleChromosome;
import io.jenetics.DoubleGene;
import io.jenetics.GaussianMutator;
import io.jenetics.Genotype;
import io.jenetics.Optimize;
import io.jenetics.Phenotype;
import io.jenetics.SinglePointCrossover;
import io.jenetics.TournamentSelector;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.engine.EvolutionStatistics;
import io.jenetics.util.ISeq;

import timeSeries.TimeSeriesDatabase;

public class Main {

	/*
	 * private static final String OUTPUT_DIRECTORY = "../harmonyresult"; private
	 * static final String TEST_NAME = "test_refactor";
	 */

	// !Achtung reihenfolge muss mit getSensorOffsetMap zusammen stimmen!
	private static final AxisStream[] axisArr = { AxisStream.GP, AxisStream.MAP, AxisStream.BP, AxisStream.SAP,
			AxisStream.WP };
	private static List<AxisStream> axisList = null;
	private static final List<String> statesToNotEvaluateList = new ArrayList<String>();

	// Test Output
	private static final boolean OUTPUT_TO_FILE = true;
	private static final String OUTPUT_DIRECTORY = "../genetic_result";
	private static final String TEST_NAME = "output_test";
	private static final int TEST_COUNT = 5;

	// The fitness function.
	private static double fitness(final Genotype<DoubleGene> genotype) {

		// Allele = Jenetic's wording for the actual gene values that
		// are stored in a chromosome (in our case of type double)
		double[] alleles = new double[genotype.geneCount()];

		int i = 0;
		for (Chromosome<DoubleGene> chromosome : genotype) {
			for (DoubleGene gene : chromosome) {
				alleles[i] = gene.doubleValue();
				i++;
			}
		}

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

		// final Factory<Genotype<DoubleGene>> gtf = Genotype.of(DoubleChromosome.of())

		final Engine<DoubleGene, Double> engine = Engine
				// Create a new builder with the given fitness
				// function and chromosome.
				.builder(Main::fitness, Genotype.of(DoubleChromosome.of(0.0, 1.0, 2), nrOfSensors))
				.offspringFraction(0.7)
				.survivorsSelector(new TournamentSelector<>(5))
				.offspringSelector(new TournamentSelector<>(5))
				.populationSize(100)
				.optimize(Optimize.MAXIMUM)
				.alterers(new GaussianMutator<>(0.1), new SinglePointCrossover<>(0.8))
				.build();

		// Exec algorithm multiple times for test
		for (int iteration = 0; iteration < TEST_COUNT; iteration++) {

			PrintStream outputFile = null;
			PrintStream resultFile = null;

			// Create output directory
			if (OUTPUT_TO_FILE) {
				String outputDir = OUTPUT_DIRECTORY + "/" + TEST_NAME;
				new File(outputDir).mkdirs();

				String fileBase = outputDir + "/" + TEST_NAME;
				try {
					outputFile = new PrintStream(
							new FileOutputStream(fileBase + "_" + iteration + "_logs" + ".csv", true), true);
					resultFile = new PrintStream(
							new FileOutputStream(fileBase + "_" + iteration + "_result" + ".txt", true), true);
					System.setOut(outputFile);
				} catch (IOException e) {
					System.err.print(e.getMessage());
					e.printStackTrace();
				}
			}

			// Create evolution statistics consumer.
			final EvolutionStatistics<Double, ?> statistics = EvolutionStatistics.ofNumber();

			StringBuilder headerStr = new StringBuilder();
			for (int i = 0; i < axisArr.length; i++) {

				headerStr.append(axisArr[i].getAxisName() + " ");
				headerStr.append(" lower;");
				headerStr.append(axisArr[i].getAxisName() + " ");
				headerStr.append(" upper;");
			}

			headerStr.append("F-Measure");

			System.out.println(headerStr);

			final ISeq<EvolutionResult<DoubleGene, Double>> best = engine.stream()
					// Truncate the evolution stream after 7 "steady"
					// generations.
					// .limit(bySteadyFitness(5))
					// The evolution will stop after maximal 1000
					// generations.
					.limit(1000)
					// Update the evaluation statistics after
					// each generation
					.peek(statistics)
					// print generation log
					.peek((evolutionResult) -> {
						Phenotype<DoubleGene, Double> bestPheno = evolutionResult.bestPhenotype();
						Genotype<DoubleGene> genum = bestPheno.genotype();
						StringBuilder str = new StringBuilder();
						for (Chromosome<DoubleGene> chromosome : genum) {
							for (DoubleGene gene : chromosome) {
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

			if (OUTPUT_TO_FILE) {
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
			if (outputFile != null) {
				outputFile.close();
			}
			if (resultFile != null) {
				outputFile.close();
			}
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

	private static void setUpDatabase(String filenameData, boolean longRun, int timespan) {
		TimeSeriesDatabase db = new TimeSeriesDatabase();
		db.setUpData(filenameData, longRun, timespan);
	}
}
