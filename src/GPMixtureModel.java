/**
 * Basic class for a univariate mixture model of a Gaussian and a Poisson distribution.
 *
 * @author tohei
 */
public class GPMixtureModel {

    /**
     * Expected value of Poisson component.
     */
    private double lambda;
    /**
     * Expected value of Gaussian component.
     */
    private double mu;
    /**
     * Standard deviation of Gaussian component.
     */
    private double sigma;

    /**
     * Mixing proportions.
     */
    private double mixProp;

    /**
     * Create new GPMixtureModel from the given parameters.
     *
     * @param lambda  expected value of the Poisson component.
     * @param mu      expected value of the Gaussian component.
     * @param sigma   standard deviation of the Gaussian component.
     * @param mixProp mixing proportions.
     */
    public GPMixtureModel(double lambda, double mu, double sigma, double mixProp) {
        this.lambda = lambda;
        this.mu = mu;
        this.sigma = sigma;
        this.mixProp = mixProp;
    }

    /**
     * Get the expected value of the Poisson component.
     *
     * @return expected value of Poisson component.
     */
    public double getLambda() {
        return lambda;
    }

    /**
     * Get the expected value of the Gaussian component.
     *
     * @return expected value of the Gaussian component.
     */
    public double getMu() {
        return mu;
    }

    /**
     * Get the standard deviation of the Gaussian component.
     *
     * @return
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Get the mixing proportions.
     *
     * @return mixing proportions.
     */
    public double getMixProp() {
        return mixProp;
    }

    /**
     * Optimize model parameters for given data vector.
     * <p>
     * This method uses the expectation maximization algorithm to optimize the set of parameters
     * for the given data.
     *
     * @param vecData SimpleVec containing data points.
     * @param eps     the algorithm terminates if the difference in log-likelihood between two iterations
     *                is below eps.
     * @param maxIter maximum number of iterations.
     * @param verbose be verbose.
     */
    public void expectationMaximization(SimpleVec vecData, double eps, int maxIter, boolean verbose) {
        int iter = 0;
        float batchSize = Math.min(100000 / vecData.length(), 1);

        if (verbose) System.out.println("Starting EM-algorithm ... ");
        while (iter < maxIter) {
            // sample from data
            SimpleVec vec = vecData.sample(batchSize);
            iter++;
            double oldLikelihood = logLikelihood(vec);
            if (verbose) {
                System.out.println("\t------------ Iteration " + iter + " ------------");
                System.out.println(toString());
                System.out.println("log-likelihood (before): " + oldLikelihood);
            }

            //expectation step
            SimpleVec y1 = normalDensity(vec, mu, sigma).times(mixProp);
            SimpleVec y2 = poissonDensity(vec, lambda).times(1 - mixProp);
            y2 = y2.add(normalDensity(vec, mu, sigma).times(mixProp));
            SimpleVec y = y1.div(y2);

            //maximization step
            SimpleVec tmp1 = (y.times(-1).add(1)); // (1-y)
            lambda = tmp1.dotProduct(vec) / MathUtils.sum(tmp1);
            mu = y.dotProduct(vec) / MathUtils.sum(y);

            SimpleVec tmp2 = vec.add(-mu);
            tmp2 = tmp2.times(tmp2);
            sigma = y.dotProduct(tmp2) / MathUtils.sum(y);

            mixProp = MathUtils.sum(y) / vec.length();

            double newLikelihood = logLikelihood(vec);
            if (verbose) System.out.println("log-likelihood (after) : " + newLikelihood);
            if (Math.abs(oldLikelihood - newLikelihood) < eps) {
                break;
            }
        }
    }

    /**
     * Compute the average log-likelihood for the current parameters under the given data.
     *
     * @param vec SimpleVec data vector.
     * @return log-likelihood of current parameters (instance variables).
     */
    public double logLikelihood(SimpleVec vec) {
        SimpleVec out = poissonDensity(vec, lambda).times(1 - mixProp);
        out = out.add(normalDensity(vec, mu, sigma).times(mixProp));
        return MathUtils.sum(MathUtils.logVec(out)) / vec.length();
    }

    /**
     * Evaluate a Poisson density function at all entries in vec.
     *
     * @param vec    SimpleVec data vector.
     * @param lambda expected value of the Poisson distribution.
     * @return new SimpleVec containing the values of the Gaussian density function evaluated at all points in vec.
     */
    public static SimpleVec poissonDensity(SimpleVec vec, double lambda) {
        double[] out = new double[vec.length()];
        for (int i = 0; i < vec.length(); i++) {
            out[i] = poissonDensity(vec.get(i), lambda);
        }
        return new SimpleVec(out);
    }

    /**
     * Evaluate a Poisson density function at the given point.
     *
     * @param x      point to evaluate.
     * @param lambda expected value of the Poisson distribution.
     * @return Poisson density at x.
     */
    public static double poissonDensity(double x, double lambda) {
        int k = (int) x;
        return Math.pow(lambda, k) * Math.exp(-lambda) / MathUtils.factorial(k);
    }

    /**
     * Evaluate a normal density function at all entries in vec.
     *
     * @param vec   SimpleVec data vector.
     * @param mu    expected value.
     * @param sigma standard deviation.
     * @return new SimpleVec containing the values of the Gaussian density function evaluated at all points in vec.
     */
    public static SimpleVec normalDensity(SimpleVec vec, double mu, double sigma) {
        double[] out = new double[vec.length()];
        for (int i = 0; i < vec.length(); i++) {
            out[i] = normalDensity(vec.get(i), mu, sigma);
        }
        return new SimpleVec(out);
    }

    /**
     * Evaluate a normal density function at a given point.
     *
     * @param x     point to evaluate.
     * @param mu    expected value.
     * @param sigma standard deviation.
     * @return normal density at x.
     */
    public static double normalDensity(double x, double mu, double sigma) {
        double c = 1 / Math.sqrt(2 * Math.PI * Math.pow(sigma, 2));
        return c * Math.exp(-Math.pow(x - mu, 2) / (2 * Math.pow(sigma, 2)));
    }

    /**
     * Generate a string representation of this mixture model.
     *
     * @return String representation of this mixture model.
     */
    @Override
    public String toString() {
        String s = "mu = " + String.format("%.4f", mu) + ", sigma = " + String.format("%.4f", sigma);
        s = s + "\nlambda = " + String.format("%.4f", lambda);
        s = s + "\nMixing proportions: " + String.format("%.3f", mixProp);
        ;
        return s;
    }
}
