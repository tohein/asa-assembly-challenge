import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

/**
 * Vector class with basic linear algebra operations.
 *
 * @author tohei
 */
public class SimpleVec {

    /**
     * Vector elements are stored in an array.
     */
    private double[] v;

    /**
     * Get the dimension of the vector.
     *
     * @return length of the vector.
     */
    public int length() {
        return v.length;
    }

    /**
     * Create a new SimpleVec from a double array.
     *
     * @param v array used to create a new SimpleVec.
     */
    public SimpleVec(double[] v) {
        this.v = v.clone();
    }

    /**
     * Create a new SimpleVec with the same data as this SimpleVec.
     *
     * @return new SimpleVec with the same entries as this object.
     */
    public SimpleVec clone() {
        return new SimpleVec(v);
    }

    /**
     * Sort this SimpleVec vector.
     */
    public void sort() {
        Arrays.sort(v);
    }

    /**
     * Compute the maximum of this SimpleVec.
     *
     * @return maximum entry of this vector.
     */
    public double max() {
        double max = v[0];
        for (int i = 0; i < this.length(); i++) {
            if (v[i] > max) max = v[i];
        }
        return max;
    }

    /**
     * Create new SimpleVec containing all entries smaller than x.
     *
     * @param x threshold.
     * @return new SimpleVec containing all entries smaller than x.
     */
    public SimpleVec smallerThan(double x) {
        LinkedList<Double> tmpList = new LinkedList<>();
        for (int i = 0; i < this.length(); i++) {
            if (v[i] < x) {
                tmpList.add(v[i]);
            }
        }
        double[] tmp = new double[tmpList.size()];
        int idx = 0;
        for (Double d : tmpList) {
            tmp[idx] = d;
            idx++;
        }
        return new SimpleVec(tmp);
    }

    /**
     * Sample with given probability from this SimpleVec.
     *
     * @param prob probability with which each element will be selected.
     * @return new SimpleVec sample from this object.
     */
    public SimpleVec sample(float prob) {
        Random rand = new Random();
        LinkedList<Double> smple = new LinkedList<>();

        for (int i = 0; i < this.length(); i++) {
            if (rand.nextFloat() < prob) {
                smple.add(v[i]);
            }
        }
        double[] tmp = new double[smple.size()];
        int idx = 0;
        for (Double d : smple) {
            tmp[idx] = d;
            idx++;
        }
        return new SimpleVec(tmp);
    }

    /**
     * Add a constant to this SimpleVec.
     *
     * @param c double constant.
     * @return new SimpleVec containing the element-wise sum of this SimpleVec and the given scalar.
     */
    public SimpleVec add(double c) {
        double[] out = new double[this.length()];
        for (int i = 0; i < this.length(); i++) {
            out[i] = this.v[i] + c;
        }
        return new SimpleVec(out);
    }

    /**
     * Add vec to this SimpleVec.
     *
     * @param vec SimpleVec to compute sum with.
     * @return new SimpleVec containing the element-wise sum of this SimpleVec and vec.
     */
    public SimpleVec add(SimpleVec vec) {
        if (this.length() != vec.length()) {
            throw new IllegalArgumentException("Incompatible dimensions.");
        } else {
            double[] out = new double[this.length()];
            for (int i = 0; i < this.length(); i++) {
                out[i] = this.v[i] + vec.v[i];
            }
            return new SimpleVec(out);
        }
    }

    /**
     * Multiply this vector with another SimpleVec.
     *
     * @param vec SimpleVec to compute product with.
     * @return element-wise product of this SimpleVec and vec.
     */
    public SimpleVec times(SimpleVec vec) {
        if (this.length() != vec.length()) {
            throw new IllegalArgumentException("Incompatible dimensions.");
        } else {
            double[] out = new double[this.length()];
            for (int i = 0; i < this.length(); i++) {
                out[i] = this.v[i] * vec.v[i];
            }
            return new SimpleVec(out);
        }
    }

    /**
     * Multiply this vector with a scalar.
     *
     * @param c a double scalar.
     * @return element-wise product of this SimpleVec and the given scalar.
     */
    public SimpleVec times(double c) {
        double[] out = new double[this.length()];
        for (int i = 0; i < this.length(); i++) {
            out[i] = this.v[i] * c;
        }
        return new SimpleVec(out);
    }

    /**
     * Divide this SimpleVec by vec.
     *
     * @param vec SimpleVec divisor.
     * @return new SimpleVec containing the element-wise division of this SimpleVec and vec.
     */
    public SimpleVec div(SimpleVec vec) {
        if (this.length() != vec.length()) {
            throw new IllegalArgumentException("Incompatible dimensions.");
        } else {
            double[] out = new double[this.length()];
            for (int i = 0; i < this.length(); i++) {
                out[i] = this.v[i] / vec.v[i];
            }
            return new SimpleVec(out);
        }
    }

    /**
     * Get the entry at position pos.
     *
     * @param pos integer index.
     * @return double scalar at position pos.
     */
    public double get(int pos) {
        return v[pos];
    }

    /**
     * Set the entry at the given index to a new scalar.
     *
     * @param pos integer index.
     * @param d   new scalar entry.
     */
    public void set(int pos, double d) {
        v[pos] = d;
    }

    /**
     * Compute dot-product between this SimpleVec and vec.
     *
     * @param vec a SimpleVec object of the same length as this vector.
     * @return the dot-product of this SimpleVec and vec.
     */
    public double dotProduct(SimpleVec vec) {
        if (this.length() != vec.length()) {
            throw new IllegalArgumentException("Incompatible dimensions.");
        } else {
            double out = 0;
            for (int i = 0; i < this.length(); i++) {
                out += this.v[i] * vec.v[i];
            }
            return out;
        }
    }

    /**
     * Creates string representation of this SimpleVec.
     *
     * @return string representation of this SimpleVec.
     */
    @Override
    public String toString() {
        return Arrays.toString(v);
    }
}
