package org.example;

import java.lang.reflect.Array;
import java.util.Arrays;


public class Main {
    public static void main(String[] args) {
        CodilityServiceMethod codilityServiceMethod = new CodilityServiceMethod();

        System.out.println(codilityServiceMethod.findingLargestBinaryGap(32));

        int[] rotateArray =  {1, 1, 2, 3, 5};
        Arrays.stream(codilityServiceMethod.rotateArrayKTime(rotateArray,2)).forEach(a -> System.out.print(a+" "));

        int[] array =  {9,3,9,3,9,7,9};

        System.out.println(codilityServiceMethod.getUnPairedValues(array));

        System.out.println(codilityServiceMethod.countNumberOfJumps(10,  85,30));

        int[] array1 =  {1,2,3,3,5,7};
        System.out.println(codilityServiceMethod.getMissingValue(array1));

        int[] intArray = {3,1,2,4,3};
        System.out.println(codilityServiceMethod.TapeEquilibrium( intArray));

        int[] intArray1 = {1, 4, 2, 3, 5, 4};
        System.out.println(codilityServiceMethod.frog(  5,intArray1));

        int[] intArray2 = {1, 4, 2,3,3, 5};
        System.out.println(codilityServiceMethod.isArrayPermutation(intArray2  ));

        int[] intArray3 = {-1,-2};
        System.out.println(codilityServiceMethod.getSmallestPositiveIntegerThatisNotPresent(intArray3));

        int[] intArray4 = {0,1,0,1,1};
        System.out.println(codilityServiceMethod.countNoOfCarPassings(intArray4));
        System.out.println(codilityServiceMethod.checkNoInRangeDivisibleByK(3,13,3));

        int[] intArray5 = {3,4,3,2,3,-1,3,3};
        System.out.println(codilityServiceMethod.dominatorOfArray(intArray5));

        int[] intArray6 = {3,2,-6,4,0};
        System.out.println(codilityServiceMethod.findMaxSumOfSubsequenceOfArrayElements(intArray6));

        String S = "CAGCCTA";
        int[] P = {2,5,0};
        int[] Q = {4,5,6};
        Arrays.stream(codilityServiceMethod.findMinimalNucleotideFromRangeOfSequenceDNA(S,P,Q)).forEach(a -> System.out.print(a+" "));

        int[] intArray7 = {4,2,2,5,1,5,8};
        System.out.println("\n"+codilityServiceMethod.findMinimalAverageOfAnySlice(intArray7));

        int[] intArray8 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.computeDistinctNumberInArray(intArray8));

        int[] intArray9 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.maximizeTriplet(intArray9));

        int[] intArray10 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.checkTriangleCanBeBuilt(intArray10));

        int[] intArray11 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.numberOfIntersectionsInDiscs(intArray11));

        String parent = "{[()()]}";
        System.out.println(codilityServiceMethod.determineParenthesesIsProper(parent));

        int[] intArray12 = {4,2,2,5,1,5,8};
        int[] intArray13 = {0,1,0,0,0};
        System.out.println(codilityServiceMethod.howManyFishAliveAlongRivers(intArray12,intArray13));

        String str = "(()(())())";
        System.out.println(codilityServiceMethod.DetermineWhetherGivenStringParentheses(str));

        int[] intArray14 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.stoneWall(intArray14));

        int[] intArray15 = {4,2,2,5,1,5,8};
        System.out.println(codilityServiceMethod.maxSlideDoubleSum(intArray15));

        System.out.println(codilityServiceMethod.countFactors(30));
        System.out.println(codilityServiceMethod.minPerimeterRectangle(30));
        System.out.println(codilityServiceMethod.flags(intArray2));

        int[] peakTestExample = {1, 2, 3, 4, 3, 4, 1, 2, 3, 4, 6, 2};
        System.out.println(codilityServiceMethod.peak(peakTestExample));
        Arrays.stream(codilityServiceMethod.countNonDivision(peakTestExample)).forEach(a -> System.out.print(a+" "));

        int[] numP =  {1, 4, 16};
        int[] numQ =  {26, 10, 20};
        Arrays.stream(codilityServiceMethod.countSemiprimes(26,numP,numQ)).forEach(a -> System.out.print(a+" "));
    }
}