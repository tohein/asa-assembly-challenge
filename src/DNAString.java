import java.util.Arrays;

/**
 * Immutable class for DNA sequences.
 *
 * @author tohei
 */
public final class DNAString implements CharSequence {

    /**
     * Valid characters.
     */
    private static final byte[] ALPHABET = {'A', 'C', 'G', 'T'};

    /**
     * Complement bases.
     */
    private static final byte[] COMPLEMENT = new byte[128];

    static {
        for (int i = 0; i < COMPLEMENT.length; i++) {
            COMPLEMENT[i] = (byte) i;
        }
        COMPLEMENT['a'] = COMPLEMENT['A'] = 'T';
        COMPLEMENT['c'] = COMPLEMENT['C'] = 'G';
        COMPLEMENT['g'] = COMPLEMENT['G'] = 'C';
        COMPLEMENT['t'] = COMPLEMENT['T'] = 'A';
    }

    /**
     * Hash value.
     */
    private int hash;

    /**
     * DNA sequence.
     */
    private final byte[] sequence;

    /**
     * Construct DNAString from byte array.
     *
     * @param sequence byte array.
     */
    public DNAString(byte[] sequence) {
        this.sequence = sequence.clone();
    }

    /**
     * Construct DNAString from String.
     *
     * @param s String sequence.
     */
    public DNAString(String s) {
        this.sequence = s.getBytes();
    }

    /**
     * Construct empty DNAString.
     */
    public DNAString() {
        this.sequence = new byte[0];
    }

    /**
     * Compute reverse complement of DNA sequence.
     *
     * @return DNAString of reverse complement.
     */
    public DNAString reverseComplement() {
        byte[] rcSequence = new byte[length()];
        for (int i = 0; i < length(); i++) {
            rcSequence[length() - i - 1] = COMPLEMENT[sequence[i]];
        }
        return new DNAString(rcSequence);
    }

    /**
     * Compute all variations at given position. If pos is outside of sequence array, return extensions.
     *
     * @param pos position index.
     * @return DNAString array containing variations of sequence at position pos.
     */
    public DNAString[] allVariations(int pos) {
        DNAString[] var = new DNAString[ALPHABET.length];
        int len = length();
        if (pos < 0) {
            pos = 0;
            len = len + 1;
        } else if (pos >= len) {
            pos = len;
            len = len + 1;
        }
        for (int i = 0; i < ALPHABET.length; i++) {
            byte[] newSeq = new byte[len];
            for (int j = 0; j < pos; j++) {
                newSeq[j] = sequence[j];
            }
            int shift = len - length();
            for (int j = pos + 1; j < len; j++) {
                newSeq[j] = sequence[j - shift];
            }
            newSeq[pos] = ALPHABET[i];
            var[i] = new DNAString(newSeq);
        }
        return var;
    }

    /**
     * Creates a new DNAString where the base at position pos will be replaced.
     * @param pos replacement index.
     * @param replacement replacement base.
     */
    public DNAString replaceCharAt(int pos, byte replacement) {
        byte[] newSequence = sequence.clone();
        newSequence[pos] = replacement;
        return new DNAString(newSequence);
    }

    /**
     * Concatenates the specified string to the end of this string.
     * <p>
     * If the length of the argument string is 0, then this DNAString object is returned. Otherwise, a new
     * DNAString object is created, representing a byte sequence that is the concatenation of the byte sequence
     * represented by this DNAString object and the byte sequence represented by the argument DNAString.
     *
     * @param s the DNAString that will will be concatenated to the end of this DNString.
     * @return a new DNAString object representing the sequence of this DNAString followed by the sequence of the
     * argument's sequence.
     */
    public DNAString concat(DNAString s) {
        if (s == null) {
            throw new NullPointerException();
        }
        if (s.length() == 0) {
            return this;
        } else {
            int len1 = this.length();
            int len2 = s.length();
            byte[] concatenation = new byte[len1 + len2];
            for (int i = 0; i < len1; i++) {
                concatenation[i] = this.sequence[i];
            }
            for (int i = len1; i < concatenation.length; i++) {
                concatenation[i] = s.sequence[i - len1];
            }
            return new DNAString(concatenation);
        }
    }

    /**
     * Returns the byte representation of the base at position pos.
     *
     * @param pos position to evaluate.
     * @return char representation of byte at position pos.
     */
    public byte byteAt(int pos) {
        return sequence[pos];
    }

    /**
     * Returns the underlying byte sequence.
     *
     * @return byte sequence of this DNAString.
     */
    public byte[] toByteArray() {
        return sequence.clone();
    }

    /**
     * Compare this DNAString to the given object. Two DNAStrings are equal if their byte sequences are the same.
     *
     * @param o an Object to compare to this DNAString.
     * @return true if the object is a DNAString object with the same byte sequence.
     */
    @Override
    public boolean equals(Object o) {
        if (o == null || !(o instanceof DNAString)) {
            return false;
        } else if (o == this) {
            return true;
        } else {
            DNAString dnaString = (DNAString) o;
            return Arrays.equals(this.sequence, dnaString.sequence);
        }
    }

    /**
     * Returns a hash code for this string. The hash code for a DNAString object is computed in the same way
     * as for a regular String object:
     * <p>
     * s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1]
     *
     * @return a hash code for this DNAString.
     */
    @Override
    public int hashCode() {
        int h = hash;
        if (h == 0 && length() > 0) {
            for (int i = 0; i < length(); i++) {
                h = 31 * h + sequence[i];
            }
            hash = h;
        }
        return h;
    }

    /**
     * Returns length of this DNAString
     *
     * @return length of the byte array containing the bases of this DNAString.
     */
    @Override
    public int length() {
        return sequence.length;
    }

    /**
     * Returns the char representation of the base at position i.
     *
     * @param i position to evaluate.
     * @return char representation of byte at position i.
     */
    @Override
    public char charAt(int i) {
        return (char) sequence[i];
    }

    /**
     * Returns a new DNAString representing a subsequence of this DNAString.
     *
     * @param from the begin index, inclusive.
     * @param to   the end index, exclusive.
     * @return the specified subsequence.
     */
    @Override
    public DNAString subSequence(int from, int to) {
        return new DNAString(Arrays.copyOfRange(sequence, from, to));
    }

    /**
     * Turn this DNAString into a regular String.
     *
     * @return a String representation of this DNAString
     */
    @Override
    public String toString() {
        return new String(sequence);
    }
}
