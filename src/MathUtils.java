
/**
 * Various mathematical functions.
 *
 * @author tohei
 */
public class MathUtils {

    /**
     * Naive implementation of the factorial function.
     *
     * @param x input integer.
     * @return x!.
     */
    public static long factorial(long x) {
        if (x < 0) throw new IllegalArgumentException("Input to factorial has to be positive.");
        if (x < 2) return 1;
        else return x * factorial(x - 1);
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
}
