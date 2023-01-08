package org.example;

import java.util.*;
import java.util.stream.IntStream;


public class CodilityServiceMethod {
    //    Lession 1
    public int findingLargestBinaryGap(int number){
        String s = Integer.toBinaryString(number);
        byte i = 0;
        int countOne = 0;
        int startOne = 0;
        int endOne = 0;
        int largestZero = 0;
        while (i< s.length())
        {
            if (s.charAt(i) == '0')
            {
                countOne++;
            }
            else if (s.charAt(i) == '1' && startOne == 0){
                startOne = 1;
            }
            else if (s.charAt(i) == '1' && startOne == 1 && endOne == 0){
                endOne = 1;
            }

            if (startOne == 1 && endOne == 1) {
                if (largestZero <= countOne) {
                    largestZero = countOne;
                }
                countOne = 0;
                endOne = 0;
            }
            i++;
        }
      return largestZero;
    }
    //    Lessson 2
    public int[] rotateArrayKTime(int[] A, int K) {

        if (K > A.length) {
            return A;
        }
        else {
            int i = 0;
            int[] newArray = new int[A.length];
            int indexToCopy = A.length - K;
            while (indexToCopy < A.length){
                newArray[i] = A[indexToCopy];
                i++;
                indexToCopy++;
            }
            indexToCopy = 0 ;
            while (indexToCopy < (A.length - K)){
                newArray[i] = A[indexToCopy];
                i++;
                indexToCopy++;
            }
            A =null;
            return newArray;
        }
    }
    public int getUnPairedValues(int[] A) {
        boolean possibleSolution = true; // to return and properly break if not possible
        for(int i = 0; i < A.length; i++) // run for all ints
        {
            possibleSolution = true; // set possible true, in case last one failed
            for(int j = 0; j < A.length; j++) // take all ints again (to compare to the rest
            {
                if (A[i] == A[j] && i != j) { // note i escape comparing to itself
                    possibleSolution = false; // if there is a math it can't be this one
                    break; // break to save resources
                }
            }
            if(possibleSolution) // if it's the solution
            {
                return A[i];
            } // return the current number (from the initial array as that is the reference number and the 2nd is for comparing)
        }
        return 0; // return default
    }
    //    Lessson 3
    public int countNumberOfJumps(int X,int Y, int D) {
        int Len= Y-X;
        if (Len%D==0)
        {
            return Len/D;
        }
        else
        {
            return (Len/D)+1;
        }
    }
    public int getMissingValue(int[] A){
        long N = A.length + 1;
        long total = N * (N + 1) / 2;
        for (int i : A)
             total -= i;
        return (int)total;
    }
    public int  TapeEquilibrium(int[] A) {
        int rightSum = 0;
        int leftSum = 0;
        int result = Integer.MAX_VALUE;

        for(int i=0; i<A.length; i++){
            rightSum += A[i];
        }

        int temp = 0;
        for(int i=1; i<A.length; i++){
            rightSum -= A[i-1];
            leftSum += A[i-1];
            temp = Math.abs(rightSum -leftSum);
            if(temp<result) result = temp;
        }

        return result;
    }
    //    Lessson 4
    public int frog(int X, int[] A) {
        int steps = X;
        boolean[] bitmap = new boolean[steps+1];
        for(int i = 0; i < A.length; i++){
            if(!bitmap[A[i]]){
                bitmap[A[i]] = true;
                steps--;
                if(steps == 0) return i;
            }

        }
        return -1;
    }
    public int isArrayPermutation(int[] A) {

        Arrays.sort(A);
        int n = A.length;
        for (int i = 1; i <= n; i++) {
            if (A[i - 1] != i) {
                return 0;
            }
        }
        return 1;
    }
    public int getSmallestPositiveIntegerThatisNotPresent(int[] A) {
        Arrays.sort(A);
        int smallest = 1;

        for (int i = 0; i < A.length; i++) {
            if (A[i] == smallest) {
                smallest++;
            }
        }
        return smallest;
    }
    //    Lessson 5
    public int countNoOfCarPassings(int[] A){
        int countOfZeros = 0, count = 0;
        for (int i = 0; i < A.length; i++){
            if (A[i] == 0) countOfZeros++;
            if (A[i] == 1) count += countOfZeros;
            if (count > 1000000000) return -1;
        }
        return count;
    }
    public int checkNoInRangeDivisibleByK(int A, int B, int K){
        int bAndKDivision = B/K;
        int a = (A > 0 ? (A - 1)/K : 0);
        if(A == 0){
            bAndKDivision++;
        }
        return bAndKDivision - a;
    }
    public int[] findMinimalNucleotideFromRangeOfSequenceDNA(String S, int[] P, int[] Q) {
        int [] C = new int[P.length];
        for (int i=0; i<P.length; i++) {
            C[i] = factor(S,P[i],Q[i]);
        }

        return C;
    }
    public int factor(String S, int i, int j) {
        if (S.substring(i,j+1).contains("A")){
            return 1;
        }else if (S.substring(i,j+1).contains("C")){
            return 2;
        }else if (S.substring(i,j+1).contains("G")){
            return 3;
        }

        return 4;
    }

    public int findMinimalAverageOfAnySlice(int[] A) {
        int minAvgIdx=0;
        double minAvgVal=(A[0]+A[1])/2;
        double currAvg;
        for(int i=0; i<A.length-2; i++){
            currAvg = ((double)(A[i] + A[i+1]))/2;
            if(currAvg < minAvgVal){
                minAvgVal = currAvg;
                minAvgIdx = i;
            }
            currAvg = ((double)(A[i] + A[i+1] + A[i+2]))/3;
            if(currAvg < minAvgVal){
                minAvgVal = currAvg;
                minAvgIdx = i;
            }
        }

        currAvg = ((double)(A[A.length-2] + A[A.length-1]))/2;
        if(currAvg < minAvgVal){
            minAvgVal = currAvg;
            minAvgIdx = A.length-2;
        }
        return minAvgIdx;
    }
    //    Lessson 6
    public int computeDistinctNumberInArray(int[] A) {
        if(A.length==0) return 0;
        if(A.length==1) return 1;
        Arrays.sort(A);
        int count=1;
        int lastElement=A[0];
        for(int i=1; i<A.length; i++)
            if(A[i]!=lastElement){
                count++;
                lastElement = A[i];
            }
        return count;
    }

    public int maximizeTriplet(int[] A) {
        Arrays.sort(A);
        int n = A.length;
        return Math.max(
                A[0] * A[1] * A[n - 1],
                A[n - 3] * A[n - 2] * A[n - 1]
        );
    }
    public int checkTriangleCanBeBuilt(int[] A){
        int n = A.length;
        if(n<3){
            return 0;
        }
        Arrays.sort(A);
        for(int i=2; i<n; i++){
            if(A[i]<(long)A[i-1]+(long)A[i-2])
                return 1;
        }
        return 0;
    }
    public int numberOfIntersectionsInDiscs(int[] A) {
        int result = 0;
        int[] dps = new int[A.length];
        int[] dpe = new int[A.length];

        for (int i = 0, t = A.length - 1; i < A.length; i++)
        {
            // The centers of discs are in the range [0, A.length - 1],
            // so if two circles intersect, then they must intersect in the range [0, A.length - 1].

            // So we can think that the start point of one disc is 0 even if this disc has negative part.
            // E.g. If there is a input: A[0] = 1
            // then this disc is start at 0-1=-1, end at 0+1=1 => [-1, 1]
            // => [0, 1]
            int s = i > A[i]? i - A[i]: 0;

            int e = t - i > A[i] ? i + A[i]: t;
            dps[s]++;
            dpe[e]++;
        }

        int activeDiscs = 0;
        for (int i = 0; i < A.length; i++)
        {
            // If there are new discs which are starting at i,
            if (dps[i] > 0)
            {
                // then there new discs are all intersecting to active discs.
                // dps[i] means there are dps[i] discs are starting
                result += activeDiscs * dps[i];

                // new started discs intersect to each other
                // E.g. if dps[i] is 3, then the count of that they intersect to each other = 3+2+1 = 3*(3-1)/2
                // Sum(1,2,...,n) = n*(n-1)/2
                result += dps[i] * (dps[i] - 1) / 2;

                // The function should return −1 if the number of intersecting pairs exceeds 10,000,000.
                if (10000000 < result) return -1;

                // of course that, these new started discs are also active discs
                activeDiscs += dps[i];
            }

            // we added new started discs, then we should have a look at that
            // if there are any disc is dying, of course, we should poll them out from active discs.
            activeDiscs -= dpe[i];
        }

        return result;
    }
    //   Lessson 7
    public int determineParenthesesIsProper(String S){
        Stack<Character> T = new Stack<Character>();
        int n = S.length();
        for (int i = 0; i < n; i++) {
            char s = S.charAt(i);
            if (s == '(' || s == '[' || s == '{') {
                T.push(s);
            } else {
                if (
                        T.size() == 0
                                || (s == ')' && T.peek() != '(')
                                || (s == ']' && T.peek() != '[')
                                || (s == '}' && T.peek() != '{')
                ) {
                    return 0;
                }
                T.pop();
            }
        }
        return T.size() == 0 ? 1 : 0;
    }
    public int howManyFishAliveAlongRivers(int[] a,int[] b) {
        int remFish = a.length;
        int i = 0;
        for (i = 0; i < b.length; i++) {
            if(b[i] != 0){
                /*remFish++; }else { */ break;
            }
        }
        Stack<Integer> myQ = new Stack<Integer>();
        for (int j = i; j < b.length; j++) {
            if(b[j] == 1)
            {
                myQ.add(j);
            }
            while(b[j] == 0 && !myQ.isEmpty()) {
                if(a[j] > a[myQ.peek()]){
                    myQ.pop(); remFish--;
                }else{
                    remFish--;
                    break;
                }
            }
        }
        return remFish;
    }
    public int DetermineWhetherGivenStringParentheses (String S) {
        LinkedList<Character> stack = new LinkedList<Character>();

        for(int i=0;i<S.length();i++) {
            char c = S.charAt(i);

            if(c == '{' || c == '[' || c == '(') {
                stack.push(c);
            } else {
                if(stack.isEmpty()) {
                    return 0;
                }

                char corresponding = stack.pop();

                if(c == ')' && corresponding != '(') {
                    return 0;
                }

                if(c == ']' && corresponding != '[') {
                    return 0;
                }

                if(c == '}' && corresponding != '{') {
                    return 0;
                }

            }
        }

        return stack.isEmpty() ? 1 : 0;
    }
    public int stoneWall(int[] H){
        Stack<Integer> stack = new Stack<>();
        int hLen = H.length;
        stack.push(H[0]);
        int count = 0;
        for (int i = 1 ; i < hLen; i++) {
            int peekVal = stack.peek();
            if (peekVal < H[i]) {
                stack.push(H[i]);
            } else if (peekVal > H[i]) {
                stack.pop();
                count++;
                if (!stack.empty()) {
                    peekVal = stack.peek();

                    while (!stack.empty() && peekVal > H[i]) {
                        count++;
                        stack.pop();
                        if (!stack.empty()) {
                            peekVal = stack.peek();
                        }
                    }
                    if (peekVal == H[i]) {
                        stack.pop();
                    }
                }
                stack.push(H[i]);
            }
        }
        count += stack.size();
        return count;
    }
    //   Lessson 8
    public int dominatorOfArray(int[] A) {
        int size = 0;
        int dominator = 0;
        for (int a : A) {
            if (size == 0) {
                size++;
                dominator = a;
            } else if (dominator != a) {
                size--;
            } else {
                size++;
            }
        }
        if (size == 0)
            return -1;

        int result = 0;
        int cnt = 0;

        for (int i = 0; i < A.length; i++) {
            if (A[i] == dominator) {
                cnt++;
                result = i;
            }
        }

        if (cnt <= A.length / 2)
            return -1;
        else
            return result;
    }
    //    Lessson 9
    public long maximumPossibleEarnings(int[] A) {
        if (A == null || A.length == 0) return 0;
        int k = 0, trackIndex = 0, profit = 0;
        int k1 = 0;
        for (int i = 1; i < A.length; i++) {
            if (k1 > k) {
                if (A[i] < A[k1]) {
                    k1 = i;
                }
                else if (A[i] - A[k1] >= profit) {
                    k = k1;
                    trackIndex = i;
                    profit = A[i] - A[k];
                }
            }
            else {
                if (A[i] >= A[trackIndex]) {
                    trackIndex = i;
                    profit = A[trackIndex] - A[k];
                }
                else if (A[i] < A[k]) {
                    k1 = i;
                }
            }
        }

        return profit;

    }
    public int findMaxSumOfSubsequenceOfArrayElements(int[] A){
        int length = A.length;
        int maxFound = A[0];
        int maxSlice = A[0];
        for(int i=1;i<length;i++) {
            maxFound = Math.max(A[i], maxFound + A[i]);
            maxSlice = Math.max(maxFound, maxSlice);
        }
        return maxSlice;
    }
    public int maxSlideDoubleSum(int A[]){
        int N = A.length;
        int[] K1 = new int[N];
        int[] K2 = new int[N];

        for(int i = 1; i < N-1; i++){
            K1[i] = Math.max(K1[i-1] + A[i], 0);
        }

        for(int i = N-2; i > 0; i--){
            K2[i] = Math.max(K2[i+1]+A[i], 0);
        }
        int max = 0;
        for(int i = 1; i < N-1; i++){
            max = Math.max(max, K1[i-1]+K2[i+1]);
        }
        return max;
    }
    //    Lessson 10
    public int countFactors(int N){
        int result = 0;
        for (int i=1; i<=(double)Math.sqrt(N); i++) {
            if(i==(double)Math.sqrt(N)) {
                result++;
            }else if(N % i == 0) {
                result = result + 2;
            }
        }
        return result;
    }
    public int minPerimeterRectangle(int N){
        int min = 1+N;
        int i=1;
        while(i*i<=N) {
            if(N % i == 0) {
                min = Math.min(min, N/i+i);
            }
            i++;
        }
        return 2*min;
    }
    public int flags(int[] A){
   ArrayList<Integer> array = new ArrayList<Integer>();
        for (int i = 1; i < A.length - 1; i++) {
            if (A[i - 1] < A[i] && A[i + 1] < A[i]) {
                array.add(i);
            }
        }
        if (array.size() == 1 || array.size() == 0) {
            return array.size();
        }
        int sf = 1;
        int ef = array.size();
        int result = 1;
        while (sf <= ef) {
            int flag = (sf + ef) / 2;
            boolean suc = false;
            int used = 0;
            int mark = array.get(0);
            for (int i = 0; i < array.size(); i++) {
                if (array.get(i) >= mark) {
                    used++;
                    mark = array.get(i) + flag;
                    if (used == flag) {
                        suc = true;
                        break;
                    }
                }
            }
            if (suc) {
                result = flag;
                sf = flag + 1;
            }else {
                ef = flag - 1;
            }
        }
        return result;

    }
    public int peak(int[] A) {
        int N = A.length;
        // Find all the peaks
        ArrayList<Integer> peaks = new ArrayList<Integer>();
        for(int i = 1; i < N-1; i++){
            if(A[i] > A[i-1] && A[i] > A[i+1]) peaks.add(i);
        }

        for(int size = 1; size <= N; size++){
            if(N % size != 0) continue; // skip if non-divisible
            int find = 0;
            int groups = N/size;
            boolean ok = true;
            // Find whether every group has a peak
            for(int peakIdx : peaks){
                if(peakIdx/size > find){
                    ok = false;
                    break;
                }
                if(peakIdx/size == find) find++;
            }
            if(find != groups) ok = false;
            if(ok) return groups;
        }
        return 0;
    }
    public int[] countNonDivision(int[] A){
        int[][] D = new int[A.length*2 + 1][2];

        for (int i = 0; i < A.length; i++) {
            D[A[i]][0]++;
            D[A[i]][1] = -1;
        }

        for (int i = 0; i < A.length; i++) {
            if (D[A[i]][1] == -1) {
                D[A[i]][1] = 0;
                for (int j = 1; j <= Math.sqrt(A[i]) ; j++) {
                    if (A[i] % j == 0 && A[i] / j != j) {
                        D[A[i]][1] += D[j][0];
                        D[A[i]][1] += D[A[i]/j][0];
                    } else if (A[i] % j == 0 && A[i] / j == j) {
                        D[A[i]][1] += D[j][0];
                    }
                }
            }
        }
        for (int i = 0; i < A.length; i++) {
            A[i] = A.length - D[A[i]][1];
        }
        return A;
    }

    public int[] countSemiprimes(int N, int[] P, int[] Q){
        Integer[] primes = sieve(N/2+1);

        int[] temp = new int[N+1];
        for (int i = 0; i < primes.length; i++) {
            for (int j = 0; j < primes.length; j++) {
                int semiPrime = primes[i] * primes[j];
                if(semiPrime <= N)
                    temp[semiPrime] = 1;
            }
        }

        int[] prefix = new int[N+1];
        for (int i = 1; i < temp.length; i++) {
            prefix[i] = temp[i] + prefix[i-1];
        }

        int[] retVal = new int[P.length];
        for (int i = 0; i < retVal.length; i++) {
            retVal[i] = prefix[Q[i]] - prefix[P[i]-1];
        }

        return retVal;
    }
    public Integer[] sieve(int n) {
        boolean[] temp = new boolean[n+1];
        for (int i = 0; i < temp.length; i++) {
            temp[i] = true;
        }
        temp[0] = temp[1] = false;

        int i = 2;
        while (i * i <= n) {
            removeProducts( temp, i );
            i++;
        }

        List<Integer> ret = new ArrayList<>();
        for (int j = 0; j < temp.length; j++) {
            if(temp[j])
                ret.add( j );
        }
        return ret.toArray( new Integer[ret.size()] );
    }
    private void removeProducts(boolean[] temp, int i) {
        for (int j = i*i; j < temp.length; j++) {
            if(temp[j] && j % i == 0) {
                temp[j] = false;
            }
        }
    }

}




