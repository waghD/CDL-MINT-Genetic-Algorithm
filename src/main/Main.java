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
import io.jenetics.DoubleGene;
import io.jenetics.Genotype;
import io.jenetics.MeanAlterer;
import io.jenetics.Mutator;
import io.jenetics.Optimize;
import io.jenetics.Phenotype;
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
	
	// The fitness function.
		private static double fitness(final double[] x) {
			double sum = 0.0;
			for(int i = 0; i < x.length; i++) {
				sum += x[i];
			}
			System.out.println(sum);
			return sum;
		}

	// Set Output directory for Test
	private static final String OUTPUT_DIRECTORY = "../harmonyresult";
	private static final String TEST_NAME = "test_refactor";

	public static void main(String[] args) {
		
		//final Factory<Genotype<DoubleGene>> gtf = Genotype.of(DoubleChromosome.of())
		
		final Engine<DoubleGene, Double> engine = Engine
				// Create a new builder with the given fitness
				// function and chromosome.
				.builder(
					Main::fitness,
					Codecs.ofVector(DoubleRange.of(0.0, 100.0)))
				.populationSize(500)
				.optimize(Optimize.MINIMUM)
				.alterers(
					new Mutator<>(0.03),
					new MeanAlterer<>(0.6))
				// Build an evolution engine with the
				// defined parameters.
				.build();

			// Create evolution statistics consumer.
			final EvolutionStatistics<Double, ?>
				statistics = EvolutionStatistics.ofNumber();

			final ISeq<EvolutionResult<DoubleGene, Double>> best = engine.stream()
				// Truncate the evolution stream after 7 "steady"
				// generations.
				.limit(bySteadyFitness(30))
				// The evolution will stop after maximal 100
				// generations.
				.limit(100)
				// Update the evaluation statistics after
				// each generation
				.peek(statistics)
				// Collect (reduce) the evolution stream to
				// its best phenotype.
				//.collect(toBestPhenotype());
				.flatMap(MinMax.toStrictlyDecreasing())

				.collect(ISeq.toISeq(10));
			System.out.println(statistics);
			//System.out.println(best.population());
			for(EvolutionResult<DoubleGene, Double> result: best) {
				System.out.println(result.population());
				System.out.println(result.population().size());
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
