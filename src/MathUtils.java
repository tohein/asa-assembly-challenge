
/**
 * Various mathematical functions.
 *
 * @author tohei
 */
public class MathUtils {

    /**
     * Factorial function based on gamma function.
     *
     * @param x input integer.
     * @return x!.
     */
    public static long factorial(long x) {
        if (x < 0) throw new IllegalArgumentException("Input to factorial has to be positive.");
        return Math.round(gamma(x + 1));
    }

    /**
     * Compute an approximation to the Gamma function at x using Lanczos formula.
     *
     * @param z input double.
     * @return gamma function at x.
     */
    static double gamma(double z) {
        // g = 7, n = 9
        z -= 1;
        final double g = 7;
        final double[] LANCZOS_COEF = {0.99999999999980993227684700473478,
                676.520368121885098567009190444019,
                -1259.13921672240287047156078755283,
                771.3234287776530788486528258894,
                -176.61502916214059906584551354,
                12.507343278686904814458936853,
                -0.13857109526572011689554707,
                9.984369578019570859563e-6,
                1.50563273514931155834e-7};

        double gam = 0.5 * Math.log(2 * Math.PI) + Math.log(z + g + 0.5) * (z + 0.5) - (z + g + 0.5);
        double A = LANCZOS_COEF[0];
        for (int i = 1; i < LANCZOS_COEF.length; i++) {
            A += LANCZOS_COEF[i] / (z + i);
        }
        return Math.exp(gam + Math.log(A));
    }

    /**
     * Sum up all the entries in the given SimpleVec object.
     *
     * @param vec SimpleVec to sum up.
     * @return sum of elements in this SimpleVec.
     */
    public static double sum(SimpleVec vec) {
        double sum = 0;
        for (int i = 0; i < vec.length(); i++) {
            sum += vec.get(i);
        }
        return sum;
    }

    /**
     * Compute the median of a SimpleVec.
     *
     * @param vec SimpleVec to compute the median of.
     * @return median of all the values in this vector.
     */
    public static double median(SimpleVec vec) {
        SimpleVec vTmp = vec.clone();
        vTmp.sort();
        if (vTmp.length() % 2 == 1) {
            return vTmp.get(vTmp.length() / 2);
        } else {
            return (vTmp.get(vTmp.length() / 2) + vTmp.get(vTmp.length() / 2 - 1)) / 2;
        }
    }

    /**
     * Compute the median absolute deviation of a SimpleVec with a scale
     * factor of 1.4826 (.75 quantile of the normal distribution).
     *
     * @param vec SimpleVec to compute MAD of.
     * @return median absolute deviation of the values in this vector.
     */
    public static double mad(SimpleVec vec) {
        return mad(vec, 1.4826);
    }

    /**
     * Compute the median absolute deviation of a SimpleVec.
     *
     * @param vec   SimpleVec to compute MAD of.
     * @param scale scale factor.
     * @return median absolute deviation of the values in this vector.
     */
    public static double mad(SimpleVec vec, double scale) {
        double m = median(vec);
        SimpleVec tmp = absVec(vec.add(-m));
        return median(tmp) * scale;
    }

    /**
     * Compute the (element-wise) absolute value of a vector.
     *
     * @param vec SimpleVec to compute absolute values of.
     * @return new SimpleVec containing the (element-wise) absolute values of a vec.
     */
    public static SimpleVec absVec(SimpleVec vec) {
        double[] out = new double[vec.length()];
        for (int i = 0; i < vec.length(); i++) {
            out[i] = Math.abs(vec.get(i));
        }
        return new SimpleVec(out);
    }

    /**
     * Compute the (element-wise) logarithm of a vector.
     *
     * @param vec SimpleVec to compute logarithm of.
     * @return new SimpleVec containing the (element-wise) logarithm of a vec.
     */
    public static SimpleVec logVec(SimpleVec vec) {
        double[] out = new double[vec.length()];
        for (int i = 0; i < vec.length(); i++) {
            out[i] = Math.log(vec.get(i));
        }
        return new SimpleVec(out);
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

    public static void main(String[] args) {
        System.out.println(factorial(5));
    }
}
