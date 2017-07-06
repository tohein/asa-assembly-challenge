
/**
 * Class for finding a suitable coverage cutoff using
 * a Gauss-Poisson mixture model.
 *
 * @author tohei
 */
public class CovCutoffFinder {

    /**
     * Gauss-Poisson mixture model.
     */
    private GPMixtureModel mixt;

    /**
     * Coverage data to analyze.
     */
    private SimpleVec covData;

    /**
     * Maximum number of tries for fitting the mixture model.
     */
    private static final int MAX_TRIES = 3;

    /**
     * Coverage cutoff will be set at the largest integer where the likelihood ratio Poisson/Gaussian > MIN_RATIO.
     */
    private static final int MIN_RATIO = 1000;

    /**
     * Create new CovCutoffFinder.
     *
     * @param covData coverage data to analyze.
     */
    public CovCutoffFinder(SimpleVec covData) {
        this.covData = covData;
    }

    /**
     * Find a suitable coverage cutoff using a Gauss-Poisson mixture model.
     *
     * @param verbose be verbose.
     * @return coverage cutoff.
     */
    public int findCutoff(boolean verbose) {
        if (verbose) System.out.println("Determine cutoff-threshold from k-mer coverage distribution: ");

        // determine median and median absolute deviation
        double median = MathUtils.median(covData);
        double mad = MathUtils.mad(covData);

        int cutoff;
        if (mad == 0) {
            // only single coverage level. no need to use cutoff.
            cutoff = 0;
        } else {
            // discard very high coverage entries (which are not very well modelled by a Gaussian and
            // therefore drive the likelihood function to zero)
            covData = covData.smallerThan(2 * median);

            // starting parameters for the mixture model
            double initVar = mad;
            double initMeanGaussian = median;
            double initMeanPoisson = 0.1;

            // use expectation maximization on a batch of size batchSize (roughly) (again, using too many
            // data points drives the likelihood function to zero)
            float batchSize = (float) 100000.0;

            // create Gaussian-Poisson mixture model
            mixt = new GPMixtureModel(initMeanPoisson, initMeanGaussian, initVar, .6);

            cutoff = 0;
            // run EM
            for (int i = 0; i < MAX_TRIES; i++) {
                float sampleProb = Math.min(batchSize / covData.length(), 1);
                mixt.expectationMaximization(covData, Math.pow(10, -4), 100, sampleProb, verbose);

                if ((mixt.getMu() < mixt.getLambda()) || (mixt.getSigma() < Math.pow(10, -4))) {
                    // improper fit, try again
                    if (verbose) System.out.println("Could not fit mixture model.");
                    batchSize = batchSize / 2;
                    mixt = new GPMixtureModel(initMeanPoisson, initMeanGaussian, initVar, .6);
                } else if (Math.abs(mixt.getMu() - mixt.getLambda()) < mixt.getSigma()) {
                    // probably single peak (no errors)
                    cutoff = 0;
                    break;
                } else {
                    cutoff = findThreshold();
                    break;
                }
            }
        }
        if (verbose) System.out.println("\nChosen cutoff threshold: " + cutoff);
        return cutoff;
    }

    /**
     * Find largest integer where the likelihood ratio Poisson/Gaussian > MIN_RATIO.
     *
     * @return estimated error threshold.
     */
    private int findThreshold() {
        int left = (int) Math.ceil(mixt.getLambda());
        for (int i = left; i >= 0; i--) {
            if (computeRatio(i) > MIN_RATIO) {
                return i;
            }
        }
        return 0;
    }

    /**
     * Compute likelihood ratio for the Gauss-Poisson mixture model.
     *
     * @param x double input value.
     * @return likelihood ratio at x.
     */
    private double computeRatio(double x) {
        double ratio = mixt.getMixProp() * MathUtils.poissonDensity(x, mixt.getLambda());
        return ratio / ((1 - mixt.getMixProp()) * MathUtils.normalDensity(x, mixt.getMu(), mixt.getSigma()));
    }
}
