
/**
 * Class for finding a suitable coverage cutoff using
 * a Gauss-Poisson mixture model.
 *
 * @author tohei
 */
public class CovCutoffFinder {

    private GPMixtureModel mixt;
    private SimpleVec covData;

    public CovCutoffFinder(SimpleVec covData) {
        this.covData = covData;
    }

    public int findCutoff(boolean verbose) {
        if (verbose) System.out.println("Determining cutoff-threshold ... ");

        double median = MathUtils.median(covData);
        double mad = MathUtils.mad(covData);
        if (mad == 0) return 0;

        // discard very high coverage entries
        covData = covData.smallerThan(2 * median);

        // create Gaussian-Poisson mixture model
        GPMixtureModel mixt = new GPMixtureModel(0.1, median, mad, .6);
        mixt.expectationMaximization(covData, Math.pow(10, -4), 100, verbose);

        return (int) Math.exp(Math.min(mixt.getLambda(), mixt.getMu()));
    }
}
