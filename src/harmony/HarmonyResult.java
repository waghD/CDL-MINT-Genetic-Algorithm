package harmony;

public class HarmonyResult {
	private boolean optimumFound;
	private int nrOfIterationsForBestFMeasure = 0;
	private int nrOfIterationsForBestMinimizedRange = 0;
	private double absOffsetBestInitial = 0.0;
	private double absOffsetBestMinimized = 0.0;
	private double runtimeAllIterations = 0;
	private double runtimeTilOptimumFound = 0;
	private double avgBestPrecision = 0;
	private double avgBestRecall = 0;
	private double avgBestFMeasure = 0;
	private double absoluteOffset = 0;
	
	public double getAbsOffsetBestInitial() {
		return absOffsetBestInitial;
	}
	
	public void setAbsOffsetBestInitial(double absOffsetBestInitial) {
		this.absOffsetBestInitial = absOffsetBestInitial;
	}
	public double getAbsOffsetBestMinimized() {
		return absOffsetBestMinimized;
	}
	
	public void setAbsOffsetBestMinimized(double absOffsetBestMinimized) {
		this.absOffsetBestMinimized = absOffsetBestMinimized;
	}
	public double getAbsoluteOffset() {
		return absoluteOffset;
	}
	
	public void setAbsoluteOffset(double offset) {
		this.absoluteOffset = offset;
	}
	public double getAvgBestPrecision() {
		return avgBestPrecision;
	}
	
	public void setAvgBestPrecision(double prec) {
		this.avgBestPrecision = prec;
	}
	
	public double getAvgBestFMeasure() {
		return avgBestFMeasure;
	}
	
	public void setAvgBestFMeasure(double fMeasure) {
		this.avgBestFMeasure = fMeasure;
	}
	
	public double getAvgBestRecall() {
		return avgBestRecall;
	}
	
	public void setAvgBestRecall(double rec) {
		this.avgBestRecall = rec;
	}
	
	public boolean getOptimumFound() {
		return optimumFound;
	}
	
	public void setOptimumFound(boolean opt) {
		this.optimumFound = opt;
	}
	public int getNrOfIterationsForBestMinimizedRange() {
		return nrOfIterationsForBestMinimizedRange;
	}
	public void setNrOfIterationsForBestMinimizedRange(int nrOfIterationsForBestMinimizedRange) {
		this.nrOfIterationsForBestMinimizedRange = nrOfIterationsForBestMinimizedRange;
	}
	public int getNrOfIterationsForBestFMeasure() {
		return nrOfIterationsForBestFMeasure;
	}
	public void setNrOfIterationsForBestFMeasure(int nrOfIterationsForBestFMeasure) {
		this.nrOfIterationsForBestFMeasure = nrOfIterationsForBestFMeasure;
	}
	public double getRuntimeIterations() {
		return runtimeAllIterations;
	}
	public void setRuntimeIterations(double runtimeIterations) {
		this.runtimeAllIterations = runtimeIterations;
	}
	
	public double getRuntimeTilOptimumFound() {
		return runtimeTilOptimumFound;
	}
	public void setRuntimeTilOptimumFound(double runtimeTilOptimumFound) {
		this.runtimeTilOptimumFound = runtimeTilOptimumFound;
	}
	
	public String toString() {
		return "Number of iterations til best f-measure found: " + this.nrOfIterationsForBestFMeasure
				+ "\nNumber of iterations til best minimized range: " + this.nrOfIterationsForBestMinimizedRange
				+ "\nAbs Offset Difference (between best f-measure and offset optimized best f-measure): " + ((1.0-(this.absOffsetBestMinimized/this.absOffsetBestInitial))*100) + "%"
				+ "\nRuntime of all iterations [s]: " + this.runtimeAllIterations
				+ "\nRuntime til optimum found [s]: " + (this.runtimeTilOptimumFound > 0.0 ? this.runtimeTilOptimumFound: "No optimum found");
		
	}
	
	
	
	
}
