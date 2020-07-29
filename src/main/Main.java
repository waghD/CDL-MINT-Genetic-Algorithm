package main;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import io.jenetics.Chromosome;
import io.jenetics.Crossover;
import io.jenetics.DoubleChromosome;
import io.jenetics.DoubleGene;
import io.jenetics.GaussianMutator;
import io.jenetics.Genotype;
import io.jenetics.MultiPointCrossover;
import io.jenetics.Optimize;
import io.jenetics.Phenotype;
import io.jenetics.SinglePointCrossover;
import io.jenetics.TournamentSelector;
import io.jenetics.UniformCrossover;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.engine.EvolutionStatistics;
import io.jenetics.stat.DoubleMomentStatistics;
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
	private static final String TEST_NAME = "states_12";
	private static final int TEST_COUNT = 10;
	private static final int populationSize = 100;


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

	public static void main(String[] args) throws FileNotFoundException {

		// Add sensors to use for detection here, derive of "AxisStream"-Class
		setUpDatabase("./lib/Daten_12_156.csv", false, 0);
		Evaluation eval = new Evaluation("./lib/realStates_12_156.csv");

		int nrOfSensors = axisArr.length;
		axisList = new ArrayList<AxisStream>(Arrays.asList(axisArr));

		eval.setUpRealDataStream(axisList);
	
		final List<Double> mutationPropabilityList = new ArrayList<Double>(Arrays.asList(0.05));
		final List<Crossover> crossoverMethodList = new ArrayList<Crossover>(
					Arrays.asList(
							new SinglePointCrossover(0.8), 
							new SinglePointCrossover(0.9),
							new MultiPointCrossover(0.8),
							new MultiPointCrossover(0.9),
							new UniformCrossover(0.8),
							new UniformCrossover(0.9)
						)
		);

		// final Factory<Genotype<DoubleGene>> gtf = Genotype.of(DoubleChromosome.of())
		for (double mutProp: mutationPropabilityList) {
			for (Crossover crossoverMethod: crossoverMethodList) {
					
						
				
				final Engine<DoubleGene, Double> engine = Engine
						// Create a new builder with the given fitness
						// function and chromosome.
						.builder(Main::fitness, Genotype.of(DoubleChromosome.of(0.0, 0.4, 2* nrOfSensors)))
						.offspringFraction(0.6)
						.survivorsSelector(new TournamentSelector<>())
						.offspringSelector(new TournamentSelector<>())
						.populationSize(populationSize)
						.optimize(Optimize.MAXIMUM)
						.alterers(new GaussianMutator<>(mutProp), crossoverMethod)
						.build();

				//Save statistics per repetition

				List<EvolutionStatistics<Double, DoubleMomentStatistics>> statisticsList = new ArrayList<EvolutionStatistics<Double, DoubleMomentStatistics>>();

				// Exec algorithm multiple times for test
				PrintStream averagedStatsFile = null;

				averagedStatsFile = new PrintStream(
						new FileOutputStream(OUTPUT_DIRECTORY + "/" + TEST_NAME + "_" + mutProp + "_"+ crossoverMethod.getClass().getSimpleName() + "_" +crossoverMethod.probability() +  "_averagedResult" + ".txt", true), true);

				long startTime = System.nanoTime();

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
									new FileOutputStream(fileBase  + "_" + mutProp + "_" + crossoverMethod.getClass().getSimpleName() + "_" +  crossoverMethod.probability() + "_" + iteration + "_logs" + ".csv", true), true);
							resultFile = new PrintStream(
									new FileOutputStream(fileBase +  "_" + mutProp + "_" + crossoverMethod.getClass().getSimpleName() + "_" + crossoverMethod.probability() + "_" + iteration + "_result" + ".txt", true), true);
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
							.limit(500)
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
							.collect(ISeq.toISeq());

					// .collect(ISeq.toISeq(10));

					if (OUTPUT_TO_FILE) {
						System.setOut(resultFile);
					}

					
					System.out.println(statistics);
					System.out.println("\n=> populationSize: " + populationSize);
					System.out.println("=> mutationPropability: " + mutProp);
					System.out.println("=> crossoverPropability: " + crossoverMethod.probability() + "\n");

			
					statisticsList.add((EvolutionStatistics<Double, DoubleMomentStatistics>) statistics);
					
					ArrayList<EvolutionResult<DoubleGene, Double>> arrayList = new ArrayList<>(best.asList());
					int reachedOptimumIteration = 0;
					double bestFitness = 0.0;
					for (int i = 0; i < arrayList.size(); i++) {
						EvolutionResult<DoubleGene, Double> curResult = arrayList.get(i);
						List<Phenotype<DoubleGene, Double>> solutions = curResult.population().asList();
						for (int j = 0; j < solutions.size(); j++) {
							Phenotype<DoubleGene, Double> curSol = solutions.get(j);
							//if (curSol.fitness() > fittestVal) {
								//System.out.println(i + ":");
								System.out.println(curSol);
								if(reachedOptimumIteration == 0 && curSol.fitness() == 1.0) {
									reachedOptimumIteration = i+1;
								}
							//}

						}
					}
					System.out.println("\n=> Found optimum in iteration " + reachedOptimumIteration);
					
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
				
				System.setOut(averagedStatsFile);
				
				System.out.println(averagedEvaluationStatistic(statisticsList));

				System.out.println("\n => Overall execution time in seconds: " +
										(System.nanoTime() - startTime) / 1000000000);
				if (averagedStatsFile != null) {
					averagedStatsFile.close();
				}
			}
		}
		
		
	}
	
	private static EvolutionStatistics<Double, DoubleMomentStatistics> averagedEvaluationStatistic(List<EvolutionStatistics<Double, DoubleMomentStatistics>> statisticsList) {
		System.out.println("+---------------------------------------------------------------------------+");
		System.out.println("|  Time statistics                                                          |");
		System.out.println("+---------------------------------------------------------------------------+");
		double selectionSum = 0, selectionMean = 0, alteringSum, alteringMean, fitnessCalculationSum, fitnessCalculationMean, overallExecutionSum, overallExecutionMean;
		double generations, altered, killed, invalids;
		EvolutionStatistics<Double, DoubleMomentStatistics> firstRes = statisticsList.get(0);
		for(int i = 1; i < statisticsList.size(); i++) {
			EvolutionStatistics<Double, DoubleMomentStatistics> curStat = statisticsList.get(i);

			firstRes.alterDuration().combine(curStat.alterDuration());
			firstRes.selectionDuration().combine(curStat.selectionDuration());
			firstRes.evolveDuration().combine(curStat.evolveDuration());
			firstRes.evaluationDuration().combine(curStat.evaluationDuration());
			firstRes.altered().combine(curStat.altered());
			firstRes.killed().combine(curStat.killed());
			firstRes.invalids().combine(curStat.invalids());
			firstRes.phenotypeAge().combine(curStat.phenotypeAge());
			firstRes.selectionDuration().combine(curStat.selectionDuration());
			firstRes.fitness().combine(curStat.fitness());

		}
		return firstRes;
	}
	
	
	private static String d(final DoubleMomentStatistics statistics) {
		return java.lang.String.format(
			"sum=%3.12f s; mean=%3.12f s",
			statistics.sum(), statistics.mean()
		);
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
