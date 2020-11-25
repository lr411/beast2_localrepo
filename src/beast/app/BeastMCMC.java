/*
* File BeastMCMC.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.app;


import beagle.BeagleFlag;
import beast.core.BEASTInterface;
import beast.core.Distribution;
import beast.core.Input;
import beast.app.beastapp.BeastDialog;
import beast.app.beastapp.BeastMain;
import beast.app.beauti.Beauti;
import beast.app.draw.ExtensionFileFilter;
import beast.app.util.Arguments;
import beast.app.util.Version;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Runnable;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.ThreadedTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.RandomTree;
//import beast.evolution.tree.RandomTree.ConstraintViolatedException;
//import //beast.evolution.tree.RandomTree.ConstraintViolatedException;
import beast.evolution.tree.Tree;
import beast.util.*;
import jam.util.IconUtils;
import javafx.util.Pair;
import beast.evolution.likelihood.*;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.json.JSONException;
import org.w3c.dom.ranges.Range;
import org.xml.sax.SAXException;

import java.lang.reflect.Field;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.filechooser.FileFilter;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.URL;
// Leo: used to set the random seed with milliseconds
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import static java.nio.file.StandardOpenOption.*;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Main application for performing MCMC runs.
 * See getUsage() for command line options.
 */
public class BeastMCMC {
    final public static String VERSION = "2.0 Release candidate";
    final public static String DEVELOPERS = "Beast 2 development team";
    final public static String COPYRIGHT = "Beast 2 development team 2011";
    // number of particles for the SMC
    public static long NR_OF_PARTICLES;
    // path to save the logs
    static String logsPath="/Users/lr411/Leo/Github/Genomics/logs_BEAST2/";
    // nr of MCMC moves
    protected long m_particleNr;
    public static final int NR_OF_MCMC_MOVES = 5;
    public BEASTInterface m_treeobj=null;
    public double m_initPopSize;
    public double m_gammaShapeLog;
    public double m_gammaShape;
    public boolean isonLeoPC=false;
    public static double targetCESSval;
    /**
     * number of threads used to run the likelihood beast.core *
     */
    static public int m_nThreads = 1;
    
    static private String[] m_sequencesArray = {
			"GTTGGCACAGTCGAATGACTGGTATACTGTTCGTCAACGATTACATAGGACTCGACTGAGCGGGACGAACTTAGGCATAATGGGGAAAGTAGCCTCCCTTCCATACCGCAAGATTTGGTATACTCTCCCCGTCTGGCAGAATCGTCCCCCTATTTTGTCCACTAGTTTACACGTGGCAGAAGCCGGACGGGGTATTGCCTCGTCTCGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAAACAAGGACCACTACAGATGTGAATGGGCACCACTTAAGCATCTGCATCGAACACATGGGAAAACCTTTCGTTCAAGCAAATCTAAACTTAGACCACCGACCTCTTGTGATGCTCATTCGGACGGAGAGTGATCCAAGGGAGTCTGGGTTACGGGTCTTGTGTAGACCTTTCAGTCAATGCCTCCATTTTGATCCAAACACTGGATGTTCAATCACTCGTTAGTAGCCCACTGTTATGAGAGGACACCAGTTGACCAGGATGCCGACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGTTGGCAAACTATCGAGACTGTTTGCACGTAAAACCGGCATTCAGCAACTGTTGACCGGTGGTAATCCTCGAGGGGTCAAATCACTGATCCATTAGAGTACCCTGTACTAGCCTATACAACCGAAGGTAAATCGATACCTAACAGGGGTTATGCGCTCCTAAACGCTTCCCAGAGTGTGCGCTGCTCGGCTAAGGCGTTCCAAACTTGTAAAAATCTTTACGCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTGTACGTTGATCGGAGCTGCGATTAACACACGGTAAGGCACGGTGGTATCAAAACTTTACATTCATGAAGTATGAGGGTGTCACAATCAAAGTACCAGGCATAAAAATTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACAAATTCGGGCGGAAGACCCTAATTTACAAAACCTCGGAGGTAAG",
            "CTGGGCTCTGGGGGTATTCTACTCGGCGCTCCGTTCGCACATCAACTGATAATAGTATCATCAGCTTACGGGTCGGCGTGGGCGGTCCATCTTGTTTGCCTAAAGTTGATAAAAGGAGGTGGCAGATGCATATTGTCGAAGGAATCAGCATGGTCTGGCATTGAAATCCAAGCCTTTATGGGCCTGTGCTGCTCCGGACCACTTGCCTCTGCGCGGATCCAGCTCGTCACTGAGACTCTCCACCAATAATTCAAGCAGCTTTGAACGTGGTAAGGATAGGACCCCGTCCCGTGGCTACGCAGTCTACTCGGACTCTATGAAGCGAAAATCAGGTTCTCATGAAGGTCTCAAATCGATGTTTCTCAATGACTATGCAATAGATGTGCGGTGTACTACGCTTCAGATGATTGAATACTTGCCCTGCATCGAGACGCAATGTGTGGTGATATGATAAAAATCAGGTCCATGGAACTTCCAGAATAGTCAACGACTGTATGCGGCGCGTCAATGAAAAACCGACTTCCGGACTAACGTGCGATTCGCCAAGGACCACACCCGAGATGGCGCACACAGAATAGACTGGGCAAATATATTATGCTACTTTTGGTCATCGGGGTACGAGAGGTAGCCTCAGAACCGGATAAGCAGCCGCCCTGCCGTAGGGGGTCTCCGTGACAATAGACTGTAATCATCACAGTCGTAATCAGGCGTTCCATACAGTTATGCTTCGCTGAGGGTCTGGTAGATGGCTTCGGACTAGACGGCAGCACCGTGTTGACGGCCTCATTACGGGTGCAAGACCGGTTTGAGCTTCACGTGTCACAGATTTTTAAGTTGCAAATCACTCATCTCCGACACAGAGGGAAGAGGTAAGCGCACTGTTCCTCTCTGACTAGACTGGAAAGGGGTTAAAGCACTTTCTCTATTGGCTCCCATATCCGTCACTTCTGACTTCATCATGATCTACACCGCAATGCCCACTTTTCTGAAAGGACCTTGG",
            "CTGGGCTAGGGGGCTATTCTACTCGGCGTTCCGTTCGCACATCAACTGATAATTGTTTCATCAGCTTACGGGTTGGTGTGGGCGGTCCATCTGGTTTGCCTAAAGGTGATAAAAGAAGTTGCCCGATGCATATTGTTGAAGCTATCATCATGGTCTGGCATTGCAATCCAAGCCTTTATGGGCCTGTGCTGCTCCGGCCCACTTGCCTCTGCGCGGATCCAGCTCGTGAATGAGACTCTCCTCCAAGAAATCAAGCAGCTCTGAACGTGGTAAGAATTGGAACACGTCCCGTAGCTACGCAGTCTACGCGGACTCTATGAAGCGAAAATCAGGTTCCCATGAAGGTCTCTAATCCGTGTTTCTCAATGAATTTGCAATAGATGTCCGGTGTACTACGCTTAAGATGATTGAAAACCTGCCCTGATTCGAGATGAAATCTGTGGTGACAAGATAGACGTTGAGTCCCTGGAACTCCCAGAATAGTCAACGACTGTATGCGGCGCGCCGATCAAAACCCGACATCCGGATTGAGGTGCGATTCGCCAAGGACCCCACCCGAGATGGGGCACACACAATATACTGGGTTAATATATTGTGCTACTTTTGGTCAACGGGGTACCAGAGGGAGCCTCTGAACCGGATAAGCAGCCCCCCCGCCGTAGGGGGTCTCCGTGACAATAGACTGTAATCATCACATTTATATTCAGGCGTACCATAGAGTTATGCTTCGCTGAGGGTCCGGTAGATGGCTTCTGTCTAGACGGCAGCAGCGTGTTGACGGCCTCATAACGGTTGCAAGACCGGTTTGAGCTCCACGGGTCACAGATTTTTAAGTTGCAAATTGCTCATCTCCGACAGAGAGGGAAGAGCTAAGGGCGCTGTCACTCTCTGACCAGACTGGAAAGGGTTTAAAGCACTGTCTCTATTGGCTCCAATCTCCGTCACTTCTGACTTCATCATGATCTACACCGCAATGCCCAGTTTTCAGAAAGGAACTTCG",
            "AAGCCCGTGCTTTACAGTCCGCATTTTACTAGTGCACTAAATACTACACTCCGTTGGAGTTCCGCCGGATAGGTATGCATCTAAAAGAGTCAGGCCTCCATTCTCTCAGCAACGCGGATTCCAGCAGAGTCTCCATTTCCGATTGTCGTCCTGATACCGGGGGAAGAGGTGTTATTTATCTCTCCACCATATTCTTGGTTGTACCCTGTCTTGGCCTTAGGGGAGCCTTACTTGCTCCTGTGACCGGGTGAACTCCCGCGTTCCCGTCCGAGGGTCTTCCTTCGAAACACGCTGTTCGTCACCACGTTCGTCCGAGAGACCCGTGCATTTCCAGACGTACTCGTGTGCTCCACGTCAATCACGGACTAATTGTAGTACAGCGCGTGTTAATTGCAGCCTTCTAAGATCCGTTAGCCAGGGGATGGATAAATCCCCTAGCGTTAGTTGAACTCAGTAGAAGACGGAACTACTACCTATATCGCTACCGACACCGGGCAGGGGGACAGATGGGGACTAGAGCCTTATATCGTAGGATGAAAAGTCCTCCCAAGACATTGTCAGACGGATCCCAGTCCATTGATAGGTGTAGCGGGTAGTCATCTAATGTGTAAGGCCAACTATGATAGTACACATCCACGCACATACGCTTCGTAGATGCCGGCCTGCCTCCTCAATCTAGTAAGGATCTACTGCATTTTATATAACAAACAACGACGGACTTGTTCCGTGCTACTATACAGTCTGGAACACAGCCATGCGTGGAAAACTCACCAGCCAGCATGGGTGTAAGGACTTCTATAGGCAGCGGTGGAAGAGTAGGGTAGTTGGTATTTCGTCAGTTGGCAAGGTATGTTAGCGGGGCGTCGAGGGTTATGACGGTTGCACAAGTCTGGGTGATTATCAACAGCAAGCGTCGTTGACCAGTACGTTCGATACCGGGAAGGTCCACGTGCGTTTACACAATGAGACTATAACCCGCGCCAAACGACACAAGAAAATA",
            "AAGACCGTGCTTTGCGGTCCACACTTAAATAGTGCAAGAATTACTACACTCCATTGGTGTTTCGCCGGATATGTATGCATCTAAAAGAGTCAGCCCTCCAGTCTCGCAGCAACCCGGATTCCAGCAGAGTCTCTATTTCGGATTGTCGTCCTGAGACCGAGGGAAGAGCTGCTCTTTATCTCTCCACCATATTCTTGGTTATACCCTGTCTTGTCCTTAGGGGAGCCTTCCTTGCGCCTGTGAACGGGAGAGCTCCCGCGTGCCCGTCCGTGGGTCTACGTTCGAAAAACTCTGTTCGTCACTACGTTCGGCCGAGAGTACCGTGCACTTCCAGACGTACTCGTGTGCTCCACGTCAGTCACTGACTAATTGTAGCACAGCGCCTTTGAATTCAAGCCTTCTAAGATCCGTTAGCCAGGGGATGGACAAATCCCCTAGCGATAGTTGAACTCAGTAGAAGACGGAACTACTACCTATATCGCTACCGACACCGGGCAGGGGGGCAGATGGGGACTACCGCCTTATATCGTAGGATGAAAAGTCCTCCCATTACAATGTCAGAGCGATCCCACCCCATTGATGGATGTAGTGGATAGTCACCTAATGTGTATGGCCAACTATGATAGTACACTTCCACGCACATACGCTTCGTTGATGCCCGCCTCCCACCTAGATCTAGTAAGTATCTCCTGCATTTTATATAACAAACAGTGACGGACTTGTTCCGTGCGACTATACAGCCTGGAAGACAGCCATGCGTGGAAAACTCACCCGCCAGCATGGGTTTAAGGACTTCTATAGGCAGCGGTGGATGAGTAGGGTAGTTGGTATTTCGTCAGTTGGCAAGGTATGTTAGCGGGGCGTCGAGGGTCATGATGGTTGCACAAGTCTCGGTGATTATCAACTCCAACCGTAGTTGACCAGTACGTTCGATACCGGGAAGGTCCAGGTGCGTTTACACAATGAGACTTTAACCCGCGCCAAACGACACAAGAAAATA",
            "AGTGCTCAAGCCGGACCTGACGCGACCAAATATCCATCTTGAGTTCCCAAGTCTCTACACACAGCGGGGAGTTCTCGCATCAACTGACCTATCGTCGCGATTATCTCAGCGGTAACCCCAGCAGTAAGAACCTAGAGATAGTCGCCGTTAAGTTGTACATTATGAGTTATTTGACAAACTTCACAAACTGCAATTCCGGCGGGCCGGACTTTCCCATTGCGCGGCTCTCTACCACGTCTGGGGAAGCACTTACATCAGTAGCTCTTGTGCTTGGCCAACACACATACGATAAAGATCCAACATCTCCGTGCGTGGGGGCAATTCCCACACACAGCATTGCATTGGTTCAGACCAGCATCTCAGAGTGCGAATAAGCGGGCAAATTTTCATTGCTACAACCGCGATCTTCGTTTATGCTCGGCCGGAAATTTGGAAAGGAGCAAAGCTGACCCACGAGCGCGAGTCCCGCTAGCAGAAAAGTCATGTTGCATGCGTAACGGCAGTACGGGCACGGGGGTCGACCGCACTACAGATGTATGCAGTAATATTTGACTAGGGCCCTCAGGTGTGTAAACAGTAAACCGGAAATTCTCTACGTTGTTTTAGTGGACTCCCCTCTCAGGTTAAGGGGCGCCGACGTAACGCGACCGGCTTTAACATTGCGATAATCAATAGGCTGCGCAATTGTAATTCTAGGTTCTAGATAAAAGTTGGATAGTGCACGTTGTAACTACCTGACTATACGCTGCAGCGTCACAAGCATAAGTCCCCTGTGGTAGTGCTCAGTAAGGCTCACTCAGGGTACGTGCAGCGTCTTTTTCGTGCAGCCGAGCATAGTCTAAACGTTTGAGTCTAAACATAGTCAGAACGGTATGCCACTTCCCTCTCGACGACTAGCCACACACCGTGTTACAGGCTGAGTCAAAAGTATTGTGCAGAAACTAAATGGCAGTACCACAAGAGTGCCTTTTTCGGGTTTACTGTGCACTTGCGAGATCGT",
            "AGTGCTCAAGACGGACGTGACGCGACCAAATATCCATCTTGAGTTCCGAAGTCTCTACACATAGCGGGGGGTTCTCGCATCAACTGACCTATCGACGCGACTATCTCTGCGGTAACCCCAGCAGTAAGAAGAAAGAGATAGTCGCCGTTAAGTTGTACATTATGAGTTATTTGACAAACTTCCCAAACTGCAATTCCGGCGGGCCGCACCTTCCCATTGCGCGGCTCTCGACCACGTCGGGGGAAGCACTTACATCAATATCTCACGTGCTTGGCAAACACACATACGATAAAGGTCCAACATCTCGGTGCGTGGGGGCAAGTCCCACACACAGCATTGCATTGGTTCAGACCAGCATCTCAGAGTGCGAATAAGCGGGCAAATTTTCTTTGCTATAACTGCGATCTTCGTGCATGCTCGGCCGGAAATTTGGAAAGGAGCAAGGCTGACCCACCAGCGCGAGTCCGGCTAGTATAAAAGTCATGTTGCATGCGTACCGTCAGTACGGGCACGGGGTTCGACCGCACAACAGGTGTTTGCAGTAATATTTGGCTAGGGCTCTCAGGTGTGTAAACAGTAAGCGGGAAATTCTCTACGTTGTTTAAGAGGACCCCCCTCTCAGGTTGAGGGTCGCCGATGTAACGCGGCCAGCTTTAACATTGCGATATTCAAGTGCCTGCGCAATCGTAATTCTAGGTTCTAGATAAAAGTTGGATAGTGCACGTTGTAACTACCTGTCTATACGCTGCAGCGTCACAAGCATAAGTCCCCAGTGGTAGTGCTCATTAAGGCTCACTCAGGGTCCGTGCAGCGTCTTTTTTGTGCAGCCGAGAATAGTTTAAACGTTTGAGTATAAACATAGTCAGAACGGTATGCCAGTCCCCTCTCGACGAGTAGCCACAGGCCGCGATAAAGGCTGAGTCAAAAGTATTGTCCAGAAACTAAATGACAGTACCACAAGAGTGCCTTTTTCCGGCGTACGGTGCAATAGCGAGATCGT",
            "GCTGGCACAGTCGAATGATTGGTATAATGTTCGTCAACCATTACATAGGACTCGTCTGAGCGGGACGAACTTAGGCATCATGGGGAAAGAAGCCTCCCTTCTATAGCGCAAGATTTGGTATACTCTCCCCGTCGGGCAGAATCGTCCCCCTATTTTTTCCACTAGTTTACAGGTGGCAGAAGCCGGACGGCGTATTGGCCCGACTAGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAGACAAGGACCACTACCGATGTGAATGGGCACCACGTAAGCATCTGAATCGCACACATGGGAAAACCTTTCGTTCAACCAAACCTAAAGTTATACCACGGGCCTTTTGTGATGCTCATTCGGACGGAGGGCCATCCAAGGGAGTCTGGGATGGGGGTCTTGGGTAGACCTTTCAGTCAAGGCCTCCATTTTGATCCAAACATCGAATGTTCAATCACTCGTTAGTAGCCGACTGTTATGAGAGGACACCAGTTGACCCGGATGCCCACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGGTGGCATCGTATCGAGACTGTTTGCACGTAAAATCGGCTTTCAGCGACTGTTGACCGGTGGTAATCCTCGAGGGGTCAATTGACTAATCTATTAGAGTACCCATTACTACTCTATATAACCGAAGGTAAATCGATGCCTAACAGGGGTTATCCGCTCCTAAACGCTTCCCAGAGTGTGCCCAGCTCGACTAAGGCGTACCAAACTTGTTAAAATCTTTACTCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTACACGTTAATCGGAGCTGCGAGTAACACACGGTAAGGCACGGTGGTATCAAACCTTTACATTCAATAAGCATTAGGGTGTCACAATGAGAGTAGCAGGCATAACATCTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACGAATTCGGGCGGAAGACACTAATTTACGACACCTTGGTGGAAAG",
            "AAGACCTTGCTTTGCGGTCCGCACTCTACTAGTGCACTAAATAGTTCACTCCGTGGGAGTGTCCCCGGATATGTATGCATCTAAAAGACTGAAAGCTTCAGTCTCTGAGCAACGCGGCTTCCAGCAGAGTCTCGTTGTCCGATTGCCGTCCTGATACCTAGGGAAGAGTTGCTCTTCCTCTCTCCACCAAATTCTTGGTTGTACCCTGTCTTTGCCTTGGGGGAGCCTTACTTGCCCCTGTGACCGGGTGAACCCACGCGTGCCCGTCCGAGGCTCTACCTACGAAAAACTCTTTTCGTCACCAGCAAGTGCCGAGGGAACCCAGCACTTCCAGACGCACTGGTGTCCACCACGTCAATACCTGACTATTTGGAGTACGGCGCCTATTAGTTTAAGACTTTGAAGATCCGCTAGCCACGGTATGGATCAATCCCCTATCCTCAGTTTAACTCAATAGAAGACGGTTCAACTCCCTATATCGCAACCGACTCCGGGCAGTGGTAGAGATGGGTACTAGCGCCTTATATCGTAGGATAAAAAGTCCTCCCATTACATTATCAGACCGTTCCCAGGCCATTGATAGATGTAGTGGGAATACACCTAATGTGTATGTCCGCCTATGATTGTACAAATCCACGCACATGCTCTGCGTGGATGCCCGCCTCCCTCCCAGATGCAGTAAGGATCTCCTGCATTTTGTTTAACAACCAACTACGGAGTTGTTCCGTGCTACACTTCAGTCTGGAAGCCAGCCATCCATCGAAAAGTCACCAGGCAGCATGGGTGGAAGGCCTTCTATAGGCAGCGGAGGAAAAGTAGGGTAGTTGGTATTTCTTTCGTTGGCAAGGTATGTTAGCGGGACGTCGAGGGCCGGGATGGTTGCACAAGTCTAGGCGATTATCAACTGCAAAAGTTGATGACCAGTCCGTTCGATACCGGGAAGGTCCCTGAGGGGTTACACATCGAGACTAGAACGCGCGCCAAACGACACAAGAAGATA",
            "AAGACCGTGCTTTGCGGTCCGCACTCTACTAGTGCACTAAATAGTTCACTCCGTGGGAGTTTCCCCGGATATGTATGCATCTAAAAGACAGAAACCGTCAGTCTCTCAGCAACGCGTCTTCCAGCAGAGTCTCGTTGTCCGATTGCCGTCCTGATACCTAGGGAAGAGTTGCTCCTCATCTCTCCACCAAATTCTTGGTTGTACCCTGTCTTTGGCTTGGGGGAGCCTTACTTGCTCCTGTGACCGGGTGAACCACCGCGTGCCCGTCCGAGGGTCTACCTTCGCAAAACTCTTTTCGTCACCAGGAACGGCAGAAGGAACCCAGCACTTCCAGACGCACTGGAGTCCTCCACGTCAATCCCTGACTATTTGTAGGACCGCGCCTATTAGTTTAAGACTTTGAAGATCCGCTAGCCACGGTATGTATAATTCCTCTATAGTTAGTTTAACTCAATAGAAGACGGTTCTACTCCCTATATCGCAACCGACTTCGGGCAGGGGAACAGATGGGTACTAGCGCCTTCTATCGTAGGATAAAAAGTCCTCCCATTACATTGTGAGGCCGTTCCCAGGCCATTGATAGATGTAGTGGGAAGACACCTAATGTGTATGATCACCTATGATTGTACAAATCCACGCACATGCTCTGCGTTGATGCCCGCCTCCCTCCCAGATCCAGTAAGGATCTCCTGCATTTTGTTTAACAACCAACTACGGAGTTGTTCCGTGCTACACTTCTGTCTGGAAGACAGCCACCCGTCGAAAAGACAACAGGCAGCATGGGTGGAAGGCCTTCTATAGGTAGCGGAGGAAAAGTAGGGTAGTTGGTATTTCTTTCGTTGGCAAGGTATGTTAGCGGGACGTCGAGGGCCGGCATGGTTGCACAAGTCTCGGCGATTATCAACTGCAAAAGTTGATGACCAGTACGTTCGATACCGGGAAGGGCCATGAGGGTTTACACATCGAGACTATAACGCGCGCCTAACGACACAAGAAGATA",
            "GTTAGCACAGTCGGATGACTGGTATAATGTTCGTCAACGATTACATAGGACTCGACTGATCGGGACGAACTTAGGCATAAGGGGGAAAGTAGCCTCCCTTCTATACCGGAAGGTTTGGTATACTCTCCCCGTCTGGCAGAATCGTCCCTCTATTTTGTCCACTAGTTTACACGTGGCAGAAGCCGGCCGGGGTATTGCCTCGACTCGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAAACAAGGACCACTACGGATGTGAGTGGGCACCACTTAAGCGTCTGCATCGAACACATGGGAAAACCTTTCGTTAAAGCAAACCTAAAGTTAGACCACGGACCTCTTGTGATGCTAATTCGGACGGAGGGCCAGCCAAGGGAGTCTGGGTTACGGGTCTTGGGTAGACCTTTCAGTCAATGCCTCCATTTTGGTCCAAACAACGAATGTTCAATCACTCGGTAGTAGCCCAGTGTTATGAGTGGACATCAGTTGACCAGGATGCCCACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGGTGGTATCCTATCGAGACTGTTTGTACGTAAAACCGGCATTCAGCAACTGTTGACCGGTGGTAATCCTCGAGGGGTCTAATCACCGATCCATTAGAGTACCCTGTACTACCCTATACAACCGAAGGTAAATCGATACCTAACAGGGGTTATGCGCTCGTAAATGCTCCCCGGAGTGTACGCTGCTCGGCTAAGGCGTTCCAAACTTGTTAAAATCTTTACTCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTATACGTTGATCGGAGCTGCGATTAACAGACGGTAACGCACGGTGGTATCAAAACATTACATTCATTAAGCATGAGGGTGTCACATTGAAAGTACCAGGCATAAAAATTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACAAATTCGGGCGGAGGACACTAATTTACAACACCTCGGAGGTAAG"
	};

    //static private boolean m_dialogInitialized=false;
    /**
     * thread pool *
     */
    public static ExecutorService g_exec = Executors.newFixedThreadPool(m_nThreads);
    /**
     * random number seed used to initialise Randomizer *
     */
    long m_nSeed = 127;
    
    /*
     * variables that are auxiliary for calculation of formulae
     * they represent position in a list
     */
    private final static int selectedLeafPosition=0;
    private final static int heightPosition=1;
    private final static int meanGaussianPosition=2;
    private final static int sigmaGaussianPosition=3;


    public static boolean adaptiveVariance;

    BeastDialog m_dialog=null;
    
    /**
     * MCMC object to execute *
     */
    Runnable m_runnable;

    protected void setParticleNr(long particleNr)
    {
    	m_particleNr=particleNr;
    }
    
    /*
     * Leo: the MCMC element that we will be using,
     * it is the same as the m_runnable member of the beastmcmc]
     * in the future the runnable will be SMC
     */
    MCMC m_mcmc;
    
    public void SetDlg(BeastDialog dialog)
    {
    	m_dialog=dialog;
    }
    
    // the following two functions initialize the randomizer (useful for random nr generation)
    public void initRandomizer(long seed)
    {    	
        Randomizer.setSeed(seed + m_particleNr<<4);
        //if(!ThreadLocalRandom.current())
        //ThreadLocalRandom.current().setSeed(seed+ (seed/2) + m_particleNr<<4);
    }
    
    // init the randomizer using the current time in milliseconds
    public void initRandomizer()
    {
	    Calendar calendar = Calendar.getInstance();
	    long randSeed = calendar.getTimeInMillis();
    	
	    initRandomizer(randSeed);
    }
    
    /**
     * parse command line arguments, and load file if specified
     * @throws IOException 
     * @throws JSONException 
     * @throws JSONParserException 
     */
    public void parseArgs(String[] args) throws IOException, XMLParserException, JSONException {
    	parseArgs(args, false);
    }    
    
    /*
     * Creates the dialog only once
     */
    static BeastDialog CreateAndShowDialog()
    {
        /*
    	final Version version = new BEASTVersion2();
        final String titleString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                "Version " + version.getVersionString() + ", " + version.getDateString() + "</p></center></html>";
        final javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");
        final String nameString = "BEAST " + version.getVersionString();
        
        final BeastDialog dialog = new BeastDialog(new JFrame(), titleString, icon);
        
       // if(m_dialogInitialized == false)
        if(true)
        {
	        if (!dialog.showDialog(nameString, 1)) {
	            return null;
	        }
        }
        
       // m_dialogInitialized=true;

        return dialog;
        */
    	
    	return null;
    }
    
    private static Arguments parseArguments(String[] args)
    {
    final Arguments arguments = new Arguments(
            new Arguments.Option[]{

//                    new Arguments.Option("verbose", "Give verbose XML parsing messages"),
//                    new Arguments.Option("warnings", "Show warning messages about BEAST XML file"),
//                    new Arguments.Option("strict", "Fail on non-conforming BEAST XML file"),
            		
            		
                    new Arguments.StringOption("exponentFile", "exponentFile", "Specify a file from where to take the exponents for SMC"),
            		new Arguments.LongOption("sleepseconds", "Specify for how many seconds you want to sleep at the beginning"),            		
            		new Arguments.Option("adaptiveVariance", "Adaptive variance of MCMC moves"),
            		new Arguments.Option("window", "Provide a console window"),
                    new Arguments.Option("options", "Display an options dialog"),
                    new Arguments.Option("working", "Change working directory to input file's directory"),
                    new Arguments.LongOption("seed", "Specify a random number generator seed"),
                    new Arguments.StringOption("prefix", "PREFIX", "Specify a prefix for all output log filenames"),
                    new Arguments.StringOption("statefile", "STATEFILE", "Specify the filename for storing/restoring the state"),
                    new Arguments.Option("overwrite", "Allow overwriting of log files"),
                    new Arguments.Option("resume", "Allow appending of log files"),
                    new Arguments.Option("validate", "Parse the XML, but do not run -- useful for debugging XML"),
                    // RRB: not sure what effect this option has
                    new Arguments.RealOption("targetCESS", "Specify the target CESS"),
                    new Arguments.IntegerOption("errors", "Specify maximum number of numerical errors before stopping"),
                    new Arguments.IntegerOption("threads", "The number of computational threads to use (default 1), -1 for number of cores"),
                    new Arguments.LongOption("nrOfSMCparticles", "The number of particles for SMC computations"),
                    new Arguments.Option("java", "Use Java only, no native implementations"),
                    new Arguments.Option("noerr", "Suppress all output to standard error"),
                    new Arguments.StringOption("loglevel", "LEVEL", "error,warning,info,debug,trace"),
                    new Arguments.IntegerOption("instances", "divide site patterns amongst number of threads (use with -threads option)"),
                    new Arguments.Option("beagle", "Use beagle library if available"),
                    new Arguments.Option("beagle_info", "BEAGLE: show information on available resources"),
                    new Arguments.StringOption("beagle_order", "order", "BEAGLE: set order of resource use"),
                    new Arguments.Option("beagle_CPU", "BEAGLE: use CPU instance"),
                    new Arguments.Option("beagle_GPU", "BEAGLE: use GPU instance if available"),
                    new Arguments.Option("beagle_SSE", "BEAGLE: use SSE extensions if available"),
                    new Arguments.Option("beagle_single", "BEAGLE: use single precision if available"),
                    new Arguments.Option("beagle_double", "BEAGLE: use double precision if available"),
                    new Arguments.StringOption("beagle_scaling", new String[]{"default", "none", "dynamic", "always"},
                            false, "BEAGLE: specify scaling scheme to use"),
                    new Arguments.Option("help", "Print this information and stop"),
                    new Arguments.Option("version", "Print version and stop"),
                    new Arguments.Option("strictversions", "Use only package versions as specified in the 'required' attribute"),
                    new Arguments.StringOption("D", "DEFINITIONS", "attribute-value pairs to be replaced in the XML, e.g., -D \"arg1=10,arg2=20\"").allowMultipleUse(),
                    new Arguments.Option("isonLeoPC", "runs on Leo pc"),
                    new Arguments.Option("sampleFromPrior", "samples from prior for MCMC analysis (by adding sampleFromPrior=\"true\" in the first run element)"),
            });

    try {
        arguments.parseArguments(args);
    } catch (Arguments.ArgumentException ae) {
    	Log.info.println();
    	Log.info.println(ae.getMessage());
    	Log.info.println();
        printUsage(arguments);
        System.exit(1);
    }
    
    return arguments;
}
    
    private static void printUsage(final Arguments arguments) {

        arguments.printUsage("beast", "[<input-file-name>]");
        Log.info.println();
        Log.info.println("  Example: beast test.xml");
        Log.info.println("  Example: beast -window test.xml");
        Log.info.println("  Example: beast -help");
        Log.info.println();
    }
    
    /**
     * parse command line arguments, and load file if specified
     * @throws IOException 
     * @throws JSONException 
     * @throws JSONParserException 
     */
  
    
    public void parseArgs(String[] args, boolean useDialog) throws IOException, XMLParserException, JSONException {
        int i = 0;
        boolean resume = false;
        boolean useStrictVersions = false;
        boolean sampleFromPrior = false;
        Map<String, String> parserDefinitions = new HashMap<>();

        File beastFile = null;
        
        Arguments arguments=parseArguments(args);

        try {
            while (i < args.length) {
                int old = i;
                if (i < args.length) {
                    if (args[i].equals("")) {
                        i += 1;
                    } else if (args[i].equals("-batch")) {
                        Logger.FILE_MODE = Logger.LogFileMode.only_new_or_exit;
                        i += 1;
                    } else if (args[i].equals("-resume")) {
                        resume = true;
                        Logger.FILE_MODE = Logger.LogFileMode.resume;
                        System.setProperty("beast.resume", "true");
                        System.setProperty("beast.debug", "false");
                        i += 1;
                    } else if (args[i].equals("-overwrite")) {
                        Logger.FILE_MODE = Logger.LogFileMode.overwrite;
                        i += 1;
                    } else if (args[i].equals("-seed")) {
                        if (args[i + 1].equals("random")) {
                            m_nSeed = Randomizer.getSeed();
                        } else {
                            m_nSeed = Long.parseLong(args[i + 1]);
                        }
                        //if(m_nSeed == 127) // Leo: default value
                        {
                        	  Calendar calendar = Calendar.getInstance();
                        	  m_nSeed=calendar.getTimeInMillis();                        	
                        }
                        i += 2;

                    } else if (args[i].equals("-threads")) {
                        m_nThreads = Integer.parseInt(args[i + 1]);
                        g_exec = Executors.newFixedThreadPool(m_nThreads);
                        i += 2;
// use BEAST environment variable to set Beast directories as colon separated list						
//					} else if (args[i].equals("-beastlib")) {
//						ClassDiscovery.setJarPath(args[i + 1]);
//						i += 2;
                    } else if (args[i].equals("-prefix")) {
                        System.setProperty("file.name.prefix", args[i + 1].trim());
                        i += 2;
                    } else if (args[i].equals("-D")) {
                        String [] strs = args[i + 1].split("=",-1);
                        for (int eqIdx = 0; eqIdx<strs.length-1; eqIdx++) {
                            int lastCommaIdx = strs[eqIdx].lastIndexOf(",");

                            if (lastCommaIdx != -1 && eqIdx == 0)
                                throw new IllegalArgumentException("Argument to -D is not well-formed: expecting comma-separated name=value pairs");

                            String name = strs[eqIdx].substring(lastCommaIdx+1);

                            lastCommaIdx = strs[eqIdx+1].lastIndexOf(",");
                            String value;
                            if (eqIdx+1 == strs.length-1) {
                                value = strs[eqIdx+1];
                            } else {
                                if (lastCommaIdx == -1)
                                    throw new IllegalArgumentException("Argument to -D is not well-formed: expecting comma-separated name=value pairs");

                                value = strs[eqIdx+1].substring(0, lastCommaIdx);
                            }
                            parserDefinitions.put(name, value);
            			}
                        i += 2;
                    } else if (args[i].equals("-strictversions")) {
                    	useStrictVersions = true;
                        i += 1;
                    } else if (args[i].equals("-sampleFromPrior")) {
                    	sampleFromPrior = true;
                        i += 1;
                    }
                    if (i == old) {
                        if (i == args.length - 1) {
                            beastFile = new File(args[i]);
                            i++;
                        } else {
                        	// Leo: taken out this
                            //throw new IllegalArgumentException("Wrong argument");
                        	i++;
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error parsing command line arguments: " + Arrays.toString(args) + "\nArguments ignored\n\n" + getUsage());
        }
        
        { // this is for debug purpose, if it runs on my PC there are local settings
        	isonLeoPC=arguments.hasOption("isonLeoPC");
        	if(isonLeoPC)
        	{
        		logsPath="/Users/lr411/Leo/Github/Genomics/logs_BEAST2/";
        	}
        	else
        	{
        		logsPath=System.getProperty("user.dir")+"/";
        	}
        	//System.out.println("Current dir is "+ logsPath);
        	
        }
        
        {
		        /*
		         * beagle info goes here
		         */
		        boolean useBeagle = arguments.hasOption("beagle") ||
		                arguments.hasOption("beagle_CPU") ||
		                arguments.hasOption("beagle_GPU") ||
		                arguments.hasOption("beagle_SSE") ||
		                arguments.hasOption("beagle_double") ||
		                arguments.hasOption("beagle_single") ||
		                arguments.hasOption("beagle_order");
		        boolean beagleShowInfo = false;
		        long beagleFlags = 0;
		        boolean useSSE = true;
		        if (arguments.hasOption("beagle_CPU")) {
		            beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
		            useSSE = false;
		        }
		        if (arguments.hasOption("beagle_GPU")) {
		            beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
		            useSSE = false;
		        }
		        if (arguments.hasOption("beagle_SSE")) {
		            beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
		            useSSE = true;
		        }
		        if (useSSE) {
		            beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
		        }
		        if (arguments.hasOption("beagle_double")) {
		            beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
		        }
		        if (arguments.hasOption("beagle_single")) {
		            beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
		        }
		
		        
		        if (beagleFlags != 0) {
		            System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));
		        }
		        if (!useBeagle) {
		            System.setProperty("java.only", "true");
		        }
        }// end of beagle part

        
        if (beastFile == null) {
            // Not resuming so get starting options...

            
           // boolean oldDialogInitialised = this.m_dialogInitialized;
            
            
          //  if(oldDialogInitialised == false)
        //    {// only do once
                // Leo: create and show only once
                //m_dialog = CreateAndShowDialog();
        	    //BeastDialog dialog=m_dialog; // just to avoid changing names

                List<String> MCMCargs = new ArrayList<>();

                if (arguments.hasOption("overwrite")) {
                    MCMCargs.add("-overwrite");
                }

                if (arguments.hasOption("resume")) {
                    MCMCargs.add("-resume");
                }
                /*
                 * Leo: taken out so that we have a version wo dialog
                switch (m_dialog.getLogginMode()) {
	                case 0: do not ovewrite 
	                    break;
	                case 1:
	                    MCMCargs.add("-overwrite");
	                    break;
	                case 2:
	                    MCMCargs.add("-resume");
	                    break;
	            }
              */
	            MCMCargs.add("-seed");
	            
	            if (arguments.hasOption("seed")) {
	            	MCMCargs.add(arguments.getLongOption("seed") + "");
	            }
	            else
	            {
	            	
	            }
	            

	            // LEO: taken dialog out
	           // MCMCargs.add(m_dialog.getSeed() + "");
	
	            if (arguments.hasOption("threads")) {
	                int threadCount = arguments.getIntegerOption("threads");
		            if (threadCount <= 0) {
		            	threadCount = Runtime.getRuntime().availableProcessors();
		            	Log.warning.println("Setting number of threads to " + threadCount);
		            }
	                MCMCargs.add("-threads");
	                MCMCargs.add(threadCount + "");
	            }
	            /*
	             * LEO: taken dialog out
	            if (m_dialog.getThreadPoolSize() > 0) {
	                MCMCargs.add("-threads");
	                MCMCargs.add(threadCount + "");
	            }
	            */
	            
          //  }
	            
            // boolean useBeagle = dialog.useBeagle();
            // LEO: taken dialog out
	            boolean useBeagle = arguments.hasOption("beagle") ||
                    arguments.hasOption("beagle_CPU") ||
                    arguments.hasOption("beagle_GPU") ||
                    arguments.hasOption("beagle_SSE") ||
                    arguments.hasOption("beagle_double") ||
                    arguments.hasOption("beagle_single") ||
                    arguments.hasOption("beagle_order");
            boolean beagleShowInfo = false;
            long beagleFlags = 0;
            /*
             * LEO: taken dialog out
            if (useBeagle) {
                beagleShowInfo = dialog.showBeagleInfo();
                if (dialog.preferBeagleCPU()) {
                    beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
                }
                if (dialog.preferBeagleSSE()) {
                    beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
                }
                if (dialog.preferBeagleGPU()) {
                    beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
                }
                if (dialog.preferBeagleDouble()) {
                    beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
                }
                if (dialog.preferBeagleSingle()) {
                    beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
                }
            }
            */

            boolean useSSE = true;
            if (arguments.hasOption("beagle_CPU")) {
                beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
                useSSE = false;
            }
            if (arguments.hasOption("beagle_GPU")) {
                beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
                useSSE = false;
            }
            if (arguments.hasOption("beagle_SSE")) {
                beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
                useSSE = true;
            }
            if (useSSE) {
                beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
            }
            if (arguments.hasOption("beagle_double")) {
                beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
            }
            if (arguments.hasOption("beagle_single")) {
                beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
            }

            
            if (beagleFlags != 0) {
                System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));
            }
            if (!useBeagle) {
                System.setProperty("java.only", "true");
            }

            File inputFile = new File(args[args.length-1]);; //dialog.getInputFile();
           
            
            if (!beagleShowInfo && inputFile == null) {
                System.err.println("No input file specified");
                System.exit(1);
            }
      //      if(oldDialogInitialised == false)
            {// only do once
            	MCMCargs.add(inputFile.getAbsolutePath());
            }
//			BeastStartDialog dlg = new BeastStartDialog();
//			if (dlg.m_bOK) {
//				parseArgs(dlg.getArgs());
//			}
            String[] argse=MCMCargs.toArray(new String[0]);
            parseArgs(argse /*MCMCargs.toArray(new String[0])*/);
            return;
        }

        Log.warning.println("File: " + beastFile.getName() + " seed: " + m_nSeed + " threads: " + m_nThreads);
        if (resume) {
            Log.info.println("Resuming from file");
        }

        if (useStrictVersions) {
        	// grab "required" attribute from beast spec
            if (beastFile.getPath().toLowerCase().endsWith(".json")) {
            	throw new IllegalArgumentException("The -strictversions flag is not implemented for JSON files yet (only XML files are supported).");
            } else {
                BufferedReader fin = new BufferedReader(new FileReader(beastFile));
                StringBuffer buf = new StringBuffer();
                String str = null;
                int lineCount = 0;
                while (fin.ready() && lineCount < 100) {
                    str = fin.readLine();
                    buf.append(str);
                    buf.append(' ');
                }
                fin.close();
                str = buf.toString();
                int start = str.indexOf("required=");
                if (start < 0) {
                	throw new IllegalArgumentException("Could not find a 'required' attribute in the XML. Add the required attribute, or run without the -strictversions flag");
                }
                char c = str.charAt(start + 9);
                start += 10;
                int end = str.indexOf(c, start);
                String packages = str.substring(start, end);
                PackageManager.loadExternalJars(packages);
            }
        } else {
            PackageManager.loadExternalJars();
        }
        

        // parse xml
        // Randomizer.setSeed(m_nSeed);
        initRandomizer();

    	
        if (beastFile.getPath().toLowerCase().endsWith(".json")) {
            m_runnable = new JSONParser(parserDefinitions).parseFile(beastFile, sampleFromPrior);
        } else {        	
        	try {
				m_runnable = new XMLParser(parserDefinitions).parseFile(beastFile, sampleFromPrior, this);
				//this.m_treeobj;
				// Leo, take here the value of the tree
			} catch (SAXException | ParserConfigurationException e) {
				throw new IllegalArgumentException(e);
			}
        }
        m_runnable.setStateFile(beastFile.getName() + ".state", resume);
    } // parseArgs

    

    public static String getUsage() {
        return "Usage: BeastMCMC [options] <Beast.xml>\n" +
                "where <Beast.xml> the name of a file specifying a Beast run\n" +
                "and the following options are allowed:\n" +
                "-resume : read state that was stored at the end of the last run from file and append log file\n" +
                "-overwrite : overwrite existing log files (if any). By default, existing files will not be overwritten.\n" +
                "-seed [<int>|random] : sets random number seed (default 127), or picks a random seed\n" +
                "-threads <int> : sets number of threads (default 1)\n" +
                "-prefix <name> : use name as prefix for all log files\n" +
                "-beastlib <path> : Colon separated list of directories. All jar files in the path are loaded. (default 'beastlib')";
    } // getUsage

    /**
     * open file dialog for prompting the user to specify an xml script file to process *
     */
    String getFileNameByDialog() {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.addChoosableFileFilter(new FileFilter() {
            @Override
			public boolean accept(File f) {
                if (f.isDirectory()) {
                    return true;
                }
                String name = f.getName().toLowerCase();
                if (name.endsWith(".xml")) {
                    return true;
                }
                return false;
            }

            // The description of this filter
            @Override
			public String getDescription() {
                return "xml files";
            }
        });

        fc.setDialogTitle("Load xml file");
        int rval = fc.showOpenDialog(null);

        if (rval == JFileChooser.APPROVE_OPTION) {
            return fc.getSelectedFile().toString();
        }
        return null;
    } // getFileNameByDialog

    public void run() throws Exception {
        g_exec = Executors.newFixedThreadPool(m_nThreads);
        m_runnable.run();
        g_exec.shutdown();
        g_exec.shutdownNow();
    } // run


    /**
     * class for starting Beast with a dialog *
     */
    class BeastStartDialog extends JDialog {
        private static final long serialVersionUID = 1L;
        boolean m_bOK = false;
        JTextField m_fileEntry;
        JTextField m_seedEntry;
        JCheckBox m_bUseGPU;
        JComboBox<String> m_mode;

        public BeastStartDialog() {
            setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
            setModalityType(DEFAULT_MODALITY_TYPE);
            init();
            setVisible(true);
        }

        String[] getArgs() {
            List<String> args = new ArrayList<>();
            args.add("-seed");
            args.add(m_seedEntry.getText());
            switch (m_mode.getSelectedIndex()) {
                case 0:
                    break;
                case 1:
                    args.add("-overwrite");
                    break;
                case 2:
                    args.add("-resume");
                    break;
            }
//			if (m_bUseGPU.isSelected()) {
//				args.add("-useGPU");
//			}
            args.add(m_fileEntry.getText());
            return args.toArray(new String[0]);
        }

        void init() {
            try {
                setTitle("Beast Start Dialog");
                Box box = Box.createVerticalBox();

                box.add(createHeader());
                box.add(Box.createVerticalStrut(10));
                box.add(createFileInput());
                box.add(Box.createVerticalStrut(10));
                box.add(Box.createVerticalBox());
                box.add(Box.createVerticalStrut(10));
                box.add(createSeedInput());
//		        box.add(Box.createVerticalStrut(10));
//		        box.add(createBeagleInput());
                box.add(Box.createVerticalStrut(10));
                box.add(createModeInput());

                box.add(Box.createVerticalGlue());
                box.add(createRunQuitButtons());
                add(box);
                int size = UIManager.getFont("Label.font").getSize();
                setSize(600 * size / 13, 500 * size / 13);
            } catch (Exception e) {
                e.printStackTrace();
                JOptionPane.showMessageDialog(this, "Could not create dialog: " + e.getMessage());
            }
        } // BeastStartDialog::init

        private Component createHeader() {
            Box box = Box.createHorizontalBox();

            String iconLocation = "beast/app/draw/icons/beast.png";
            ImageIcon icon = null;
            try {
                URL url = ClassLoader.getSystemResource(iconLocation);
                if (url == null) {
                    System.err.println("Cannot find icon " + iconLocation);
                    return null;
                }
                icon = new ImageIcon(url);
            } catch (Exception e) {
                System.err.println("Cannot load icon " + iconLocation + " " + e.getMessage());
                return null;
            }


            JLabel label = new JLabel(icon);
            label.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 10));
            box.add(label, BorderLayout.WEST);
            label = new JLabel("<html><center>BEAST<br>Version: " + VERSION + "<br>Developers: " + DEVELOPERS + "<br>Copyright: " + COPYRIGHT + "</html>");
            label.setHorizontalAlignment(SwingConstants.CENTER);
            box.add(label);
            return box;
        } // BeastStartDialog::createHeader

        private Component createFileInput() {
            Box box = Box.createHorizontalBox();
            box.add(new JLabel("Beast XML File: "));
            m_fileEntry = new JTextField();
            int fontsize = m_fileEntry.getFont().getSize();
            Dimension size = new Dimension(300 * fontsize / 13, 20 * fontsize / 13);
            m_fileEntry.setMinimumSize(size);
            m_fileEntry.setPreferredSize(size);
            m_fileEntry.setSize(size);
            m_fileEntry.setToolTipText("Enter file name of Beast 2 XML file");
            m_fileEntry.setMaximumSize(new Dimension(1024 * fontsize / 13, 20 * fontsize / 13));
            box.add(m_fileEntry);
            //box.add(Box.createHorizontalGlue());

            JButton button = new JButton("Choose file");
            button.addActionListener(e -> {
                    JFileChooser fileChooser = new JFileChooser(Beauti.g_sDir);
                    File file = new File(m_fileEntry.getText());
                    if (file.exists())
                        fileChooser.setSelectedFile(file);
                    fileChooser.addChoosableFileFilter(new ExtensionFileFilter(".xml", "Beast xml file (*.xml)"));
                    fileChooser.setDialogTitle("Select Beast 2 XML file");
                    int rval = fileChooser.showOpenDialog(null);
                    if (rval == JFileChooser.APPROVE_OPTION) {
                        String fileName = fileChooser.getSelectedFile().toString();
                        if (fileName.lastIndexOf('/') > 0) {
                            Beauti.g_sDir = fileName.substring(0, fileName.lastIndexOf('/'));
                        }
                        m_fileEntry.setText(fileName);
                    }
                });
            box.add(button);

            return box;
        } // BeastStartDialog::createFileInput

        private Component createSeedInput() {
            Box box = Box.createHorizontalBox();
            box.add(new JLabel("Random number seed: "));
            m_seedEntry = new JTextField("127");
            m_seedEntry.setHorizontalAlignment(SwingConstants.RIGHT);
            Dimension size = new Dimension(100, 20);
            m_seedEntry.setMinimumSize(size);
            m_seedEntry.setPreferredSize(size);
            m_seedEntry.setSize(size);
            m_seedEntry.setToolTipText("Enter seed number used for initialising the random number generator");
            int fontsize = m_seedEntry.getFont().getSize();
            m_seedEntry.setMaximumSize(new Dimension(1024 * fontsize / 13, 20 * fontsize / 13));
            box.add(m_seedEntry);
            box.add(Box.createHorizontalGlue());
            return box;
        } // BeastStartDialog::createSeedInput

        private Component createModeInput() {
            Box box = Box.createHorizontalBox();
            box.add(new JLabel("Mode of running: "));
            m_mode = new JComboBox<>(new String[]{"default: only write new log files",
                    "overwrite: overwrite log files",
                    "resume: appends log to existing files (if any)"});
            Dimension size = new Dimension(350, 20);
            m_mode.setMinimumSize(size);
            m_mode.setPreferredSize(size);
            m_mode.setSize(size);
            m_mode.setMaximumSize(size);

            m_mode.setSelectedIndex(0);
            box.add(m_mode);
            box.add(Box.createHorizontalGlue());
            return box;
        } // BeastStartDialog::createModeInput

        Component createRunQuitButtons() {
            Box cancelOkBox = Box.createHorizontalBox();
            cancelOkBox.setBorder(new EtchedBorder());
            JButton okButton = new JButton("Run");
            okButton.addActionListener(e -> {
                    m_bOK = true;
                    dispose();
                });
            JButton cancelButton = new JButton("Quit");
            cancelButton.addActionListener(e -> {
                    dispose();
                    System.exit(0);
                });
            cancelOkBox.add(Box.createHorizontalGlue());
            cancelOkBox.add(cancelButton);
            cancelOkBox.add(Box.createHorizontalStrut(20));
            cancelOkBox.add(okButton);
            cancelOkBox.add(Box.createHorizontalStrut(20));
            return cancelOkBox;
        } // BeastStartDialog::createRunQuitButtons

    } // class BeastStartDialog

    // this get max val is obsolete
    public static double getMaxValue(double[] numbers){
    	  double maxValue = numbers[0];
    	  for(int i=1;i < numbers.length;i++){
    	    if(numbers[i] > maxValue){
    		  maxValue = numbers[i];
    		}
    	  }
    	  return maxValue;
    	}

    public static double logSumOfExponentials(double[] xs, double max) {
        //if (xs.length == 1) return xs[0];
        //double max = maximum(xs);
        double sum = 0.0;
        for (int i = 0; i < xs.length; ++i)
            if (xs[i] != Double.NEGATIVE_INFINITY)
                sum += java.lang.Math.exp(xs[i] - max);
        return max + java.lang.Math.log1p(sum);
    }

    public static double logSumOfExponentials(double[] xs) {
        OptionalDouble max = Arrays.stream(xs).max();
        
        return logSumOfExponentials(xs, max.getAsDouble());
    }

    // this function takes as input the normalized weights (log format)
    // and provides as output the indices of the particles that have survived the resample
    public static List<Integer> stratified_resample(double[] log_weights)
	{
      // the length of the weights is the number of particles
      int P = log_weights.length, N=P;
      double N_double=(double)N;
      // following list has fixed dimension
      List<Integer> resampleList = Arrays.asList(new Integer[N]);

      // define and set the array of weights (exponential of log weights)
      double weights[] = new double[P]; // Arrays.copyOf(log_weights, P);
      Arrays.parallelSetAll(weights, e->java.lang.Math.exp(log_weights[e]));
   
      // create variables for cumulative sum
      double sumW=Arrays.stream(weights).sum();
      double cw[] = Arrays.copyOf(weights, P);
      Arrays.parallelPrefix(cw, Double::sum);
      Arrays.parallelSetAll(cw, e-> cw[e]/sumW);
      
      Integer[] range = IntStream.rangeClosed(0, N-1).boxed().collect(Collectors.toSet()).toArray(new Integer[N]);
     
      double v[] = new double[N];
      Arrays.parallelSetAll(v, i -> (range[i] + (ThreadLocalRandom.current().nextDouble()))/N_double);

	      {
    	      int i=0,j=0;
    	      
    	      for(;i<N;i++)
    	      {
    	    	  while(cw[j]<v[i])
    	    	  {
    	    		  j++;
    	    	  }
    	    	  resampleList.set(i, j);
    	      }
	      }
      
      return resampleList;
	}// end of stratified resample

/*
score <- function(current_log_weights,target_at_particles,IS_proposal_at_particles,current_temp,previous_temp,alpha)
{
  N = length(current_log_weights)
  incremental_log_weights = (target_at_particles)*(current_temp-previous_temp)+(IS_proposal_at_particles)*(previous_temp-current_temp)
  CESS(current_log_weights,incremental_log_weights)-alpha*N
}

bisection <- function(range,current_log_weights,target_at_particles,IS_proposal_at_particles,alpha)
{
  current_bisect_temp = range[1]
  bisect_size = (range[2] - range[1])/2
  direction = 1
  
  old_score = score(current_log_weights,target_at_particles,IS_proposal_at_particles,current_bisect_temp,current_bisect_temp,alpha)
  
  for (i in 1:100)
  {
    new_bisect_temp = current_bisect_temp + direction*bisect_size
    bisect_size = bisect_size/2
    the_score = score(current_log_weights,target_at_particles,IS_proposal_at_particles,new_bisect_temp,range[1],alpha)
    direction = sign(the_score)
    
    current_bisect_temp = new_bisect_temp
    
    if (the_score==0)
    {
      break
    }
    
  }
  
  current_bisect_temp
}

IS_ESS = function(log_weights)
{
  exp(-logsumexp(2*log_weights))
}

 incremental_log_weights = (target_at_particles)*(1-previous_temp)+(IS_proposal_at_particles)*(previous_temp-1)
 log_weights = log_weights + incremental_log_weights
 log_sum_weights = logsumexp(log_weights)
 log_marginal_likelihood = log_marginal_likelihood + log_sum_weights
 log_weights = log_weights - log_sum_weights # Normalise the weig

    effective_number_of_particles = IS_ESS(log_weights)

	CESS = function(current_log_weights,incremental_log_weights)
			{
			  N = length(current_log_weights)
			  exp(log(N) + 2*logsumexp(current_log_weights+incremental_log_weights) - logsumexp(current_log_weights+2*incremental_log_weights))
			}

*/
   
    // the following are normalised and unnormalised versions of ESS calculations
    public static double CESS_Normalised(double[] incremental_log_weights, double[] normalized_log_weights)
    {
    	int N= normalized_log_weights.length;
    	
    	return CESS(incremental_log_weights, normalized_log_weights)/N;
    }

    public static double CESS(double[] incremental_log_weights, double[] normalized_log_weights)
    {
    	int N= normalized_log_weights.length;
    	// auxiliary vector used for calc
    	double[] auxiliaryVect= new double[N];
    	
    	// do incremental_log_weights + normalized_log_weights
    	Arrays.parallelSetAll(auxiliaryVect, e->incremental_log_weights[e]+normalized_log_weights[e]);
    	// numerator in log is log(N) + 2*logsumexp(normalized_log_weights+incremental_log_weights)
    	double logNumerator=java.lang.Math.log(N)+2*logSumOfExponentials(auxiliaryVect);
    	
    	// do 2*incremental_log_weights + normalized_log_weights
    	Arrays.parallelSetAll(auxiliaryVect, e->(2*incremental_log_weights[e])+normalized_log_weights[e]);
    	double logDenominator=logSumOfExponentials(auxiliaryVect);
    	
    	return Math.exp(logNumerator-logDenominator);
    }

    // the following are normalised and unnormalised versions of ESS calculations
    public static double ESS_Normalised(double[] normalized_log_weights)
    {
    	int N= normalized_log_weights.length;
    	
    	return ESS(normalized_log_weights)/N;
    }
    
    public static double ESS(double[] normalized_log_weights)
    {
    	int N= normalized_log_weights.length;
    	// auxiliary vector used for calc
    	double[] weightsSquared= new double[N];//Arrays.copyOf(normalized_log_weights, normalized_log_weights.length);
    	// multiply by 2
    	Arrays.parallelSetAll(weightsSquared, e->2.0*normalized_log_weights[e]);
    	
    	return Math.exp(-logSumOfExponentials(weightsSquared));
    }
    
/*   
     getTreePosition looks for the position of the tree
    in a state vector array
    return values:
    	position (>=0) if found
    	-1 if not found
*/
    public static int getTreePosition(MCMC mcmc)
    {
     	   int treeposition;
    	
    	   State curstate=mcmc.getState();
     	   int statenodeslength = curstate.stateNode.length;
     	   StateNode stn=null;
	         int i;  
     	     for(i=0; i< statenodeslength;i++)
	           	{
	           		stn=curstate.stateNode[i];
	           		//Class cls=stn.getClass();
	           		//if(cls.getName().endsWith(".Tree"))
	           		if(stn instanceof Tree)
	           		{
	           			treeposition=i;
	           			break;
	           		}
	           	}
	           
		        if(i<statenodeslength)
	            {
	            	treeposition=i;
	            }
	            else
	            {
	            	treeposition=-1;
	            }
		return treeposition;
    }

    /*   
    getPopsizePosition looks for the position of the parameter estimating the effective population size
   in a state vector array
   return values:
   	position (>=0) if found
   	-1 if not found
*/
   public static int getPopsizePosition(MCMC mcmc)
   {
 	   int position;
	
	   State curstate=mcmc.getState();
 	   int statenodeslength = curstate.stateNode.length;
 	   StateNode stn=null;
         int i;  
 	     for(i=0; i< statenodeslength;i++)
           	{
           		stn=curstate.stateNode[i];
           		if(stn.getID().contains("popSize"))
           		{
           			position=i;
           			break;
           		}
           	}
           
	        if(i<statenodeslength)
            {
            	position=i;
            }
            else
            {
            	position=-1;
            }
	return position;    	
}
   
/*    
 * deep copies two states from two Sequential (or BeastMCMC elements)
 * it is assumed that both elements are not null
*/    
    public static void copyState(StateNode[] source, Sequential sink)
    {    	
    	// int statenodeslength = source.m_mcmc.getState().stateNode.length;
    	//StateNode[] sourceStatenode=source.m_mcmc.getState().stateNode;
    	StateNode[] sourceStatenode=source;
    	StateNode[] sinkStatenode=sink.m_mcmc.getState().stateNode;
 
        Arrays.parallelSetAll(sinkStatenode, e ->
	       	{ 
	       		sinkStatenode[e].assignFromFragile(sourceStatenode[e]);
	       		return sinkStatenode[e];
	       	}
       	);
    }
    
    public static void copyState(Sequential source, Sequential sink)
    {    	
    	// int statenodeslength = source.m_mcmc.getState().stateNode.length;
    	StateNode[] sourceStatenode=source.m_mcmc.getState().stateNode;
    	//StateNode[] sourceStatenode=source;
    	StateNode[] sinkStatenode=sink.m_mcmc.getState().stateNode;
 
        Arrays.parallelSetAll(sinkStatenode, e ->
	       	{ 
	       		if(e==0)
	       		{
	       			int ellade;
	       			ellade=3;
	       		}
	       		sinkStatenode[e].assignFromFragile(sourceStatenode[e]);
	       		return sinkStatenode[e];
	       	}
       	);
    }

    public static int getStateSpaceDimension(Sequential source)
    {
    	int statenodeslength = source.m_mcmc.getState().stateNode.length;
    	return statenodeslength;
    }
    
	public static void updateParticlesList(List<Integer> stratifiedList, Sequential[] beastMClist)
	{
		int localIndex;
		int N_int=beastMClist.length;
	
		// process the resample: set particles according to the stratified list of particles who made it
	    for(int i=N_int-1;i>=0;i--)
		{
	    	localIndex = stratifiedList.get(i);
	
	    	if(localIndex!=i)
	    	{// only if the former position is different from the current
	    		copyState(beastMClist[localIndex], beastMClist[i]);	    		
	    	}
		}
	
	}

	// does the mcmc move and also set the exponent for the following annealing
	// [Leo: consider to split these two moves in the future for better readability of the code
	// at the moment we keep them together for efficiency]
    public static void doMCMC_andSetExponentForAnnealing(Sequential[] beastMClist, final double currentExponent, long[] nrOfMCMCrejections)
    {
       	Arrays.parallelSetAll(beastMClist, e->{
			//BeastMCMC bmc=beastMClist[e];
       		Sequential bmcc=beastMClist[e];
       	    // here update the weights
        	Tree tretre=(Tree)beastMClist[e].m_mcmc.getState().stateNode[0];

		//	System.out.println("Before particle "+e+", nodes: "+tretre.getNodeCount()+", stored: "+tretre.getStoredNodes().length);
			
       	    MCMC mc=bmcc.m_mcmc;

       	    // the mcmc run is done with the previous exponent
        	try {
        		mc.setInitState(false);
				bmcc.run();
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

		//	System.out.println("After particle "+e+", nodes: "+tretre.getNodeCount()+", stored: "+tretre.getStoredNodes().length);
        	// here get the nr of rejections
        	nrOfMCMCrejections[e]=mc.getNrOfMCMCrejections();
        	// setting the exponent sets the target distribution
        	mc.setSimulatedAnnhealingExponent(currentExponent);
       		
       		return bmcc;
			});
    }

    // initialises the normalised weights to 1/N
    public static void initNotmalisedWeights(double[] logNormalisedWeights, double minuslogN)
    {
       	// in the following initialize the weights to 1/N
		Arrays.parallelSetAll(logNormalisedWeights, e->{return minuslogN;});
    }
    
    // the output is the updated incremental weights after annealing
    public static void calculateIncrementalWeightsForTransformation(Sequential[] beastMClist, double[] output_logUnnormalisedIncrementalWeights, final double annealingExponent, boolean calledBeforeTreeUpdated)
    {
       	Arrays.parallelSetAll(output_logUnnormalisedIncrementalWeights, e->{
			//BeastMCMC bmc=beastMClist[e];
			MCMC mc=beastMClist[e].m_mcmc;
			// performance wise, we don't need the log1 and log 2 we could substitute the full expressions to log1 and 2
			double log=mc.calculateLogPSimulatedAnnhealing(annealingExponent);
//        	// ?? is it ok to set the exponent here????
			if(Double.isNaN(log))
			{
	            throw new RuntimeException(
	                    "Bad value in log expression of calculateIncrementalWeightsForTransformation\n");
			}
			
			// in the first call it calculates the old pi with old tree (before sequence add), in second call calculates the new pi with new tree (with sequence added)
			// so if it is the second call, the vector output_logUnnormalisedIncrementalWeights already contains the calculation of old pi with old tree (before sequence add)
			return calledBeforeTreeUpdated?(log):log-output_logUnnormalisedIncrementalWeights[e];
			});
    }

    public static void calculateIncrementalWeightsForTransformationAddLeafComponent(Sequential[] beastMClist, double[] output_logUnnormalisedIncrementalWeights, ArrayList<List<Double>> values)
    {
/*       	Arrays.parallelSetAll(output_logUnnormalisedIncrementalWeights, e->{
       	    return output_logUnnormalisedIncrementalWeights[e]-values.get(e).get(logLeafComponentPosition);
       	});
*/    }
    
    public static void calculateIncrementalWeightsForTransformationAddHeightComponent(Sequential[] beastMClist, double[] output_logUnnormalisedIncrementalWeights, ArrayList<List<Double>> values)
    {
/*       	Arrays.parallelSetAll(output_logUnnormalisedIncrementalWeights, e->{
       	    return output_logUnnormalisedIncrementalWeights[e]-values.get(e).get(logHeightComponentPosition);
       	});
*/    }

    // the output is the updated incremental weights after annealing
    public static void calculateIncrementalWeights(Sequential[] beastMClist, double[] output_logUnnormalisedIncrementalWeights, final double previousExponent, final double nextExponent)
    {
       	Arrays.parallelSetAll(output_logUnnormalisedIncrementalWeights, e->{
			//BeastMCMC bmc=beastMClist[e];
			MCMC mc=beastMClist[e].m_mcmc;
			// performance wise, we don't need the log1 and log 2 we could substitute the full expressions to log1 and 2
			double log1=mc.calculateLogPSimulatedAnnhealing(nextExponent);
			double log2=mc.calculateLogPSimulatedAnnhealing(previousExponent);
//        	// ?? is it ok to set the exponent here????
			if(Double.isNaN(log1-log2))
			{
	            throw new RuntimeException(
	                    "Bad value in log expression of calculateIncrementalWeights\n");
			}
			return (log1-log2);
			});
    }
    
    public static void normaliseWeights(double[] inputLogIncrementalWeights, double[] outputLogNormalizedWeights)
    {
        // multiply (i.e. add in log space) incremental part by the normalised weights to have updated unnormalised weights
		// temporarily we have the outputLogNormalizedWeights containing the multiplication (i.e. log-sum)
    	Arrays.parallelSetAll(outputLogNormalizedWeights, e->inputLogIncrementalWeights[e]+outputLogNormalizedWeights[e]);
    	
    	// normalising constant below
       	final double normalizingConstant=logSumOfExponentials(outputLogNormalizedWeights);
       
        // we normalize the weights
		Arrays.parallelSetAll(outputLogNormalizedWeights, e->outputLogNormalizedWeights[e]-normalizingConstant);
    }
    
    public static void initParticlesAndSampleFromPrior(Sequential[] beastMClist, BeastDialog dlg, String[] args)
    {
        Arrays.parallelSetAll(beastMClist, e ->
       	{ // here probably better not to call the deepcopy method we created
       		// because we need to sample from prior
       		Sequential bmc= new Sequential();
            //if(!ThreadLocalRandom.current())
            //  ThreadLocalRandom.current().setSeed(e);

       		bmc.setParticleNr(e);
            //bmc.SetDlg(dlg);
       	    // the mcmc run is done with the previous exponent
        	try {
                bmc.parseArgs(args);
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
            MCMC mc= (MCMC)bmc.m_runnable;
            // set particle nr also in the mcmc
            mc.setParticleNr(e);
        	RealParameter popParam=(RealParameter) mc.getState().stateNode[2]; //.values[0]=newParamVal;
        	StateNode[] stnl=mc.getState().stateNode;
        	double paramthis=popParam.getValue();
            mc.SetDistributionsFromInput();
            // initialise the state of the posterior
            mc.initStateAndPosterior();
            //aaaaaaaaaaaaaaaaa
        	double paramthisAfter=popParam.getValue();
            
            // this is to avoid repeated casting afterwards
            bmc.m_mcmc=(MCMC) bmc.m_runnable;
            
            return bmc;
       	});
    }
    
    public static String getDateString()
    {
		Date date = Calendar.getInstance().getTime();  
	    DateFormat dateFormat = new SimpleDateFormat("yyyymmdd_hhmmss");  
	    String strDate = dateFormat.format(date);
	    return strDate;
    }
    
    public static String formatAppendString(int N_int, long executionTime, int nrOfSteps)
    {
	    String appendString="_"+getDateString()+"_P"+N_int+"_E"+executionTime+"_S"+nrOfSteps;
	    
	    return appendString;
    }
    
    public static long getExecutionLength(long timeInNanoseconds)
    {
    	//final divider=60000000000;//1000000000*60;// 10^6 is milliseconds, 10^9 seconds, 60*10^9 minutes
    	// ms
    	long executionTime=timeInNanoseconds/1000000;
    	// s
    	executionTime=executionTime/1000; 
    	// m
    	//executionTime=executionTime/60; 

    	return executionTime;
    }
    
    public static String createAppendString(String filename,String appendString ,String fileExt)
    {
    	String outString=logsPath+filename+appendString+"."+fileExt;
    	return outString;
    }
    
    public static String createTxtAppendString(String filename,String appendString )
    {
    	return createAppendString(filename, appendString ,"txt");
    }
    
    public static String createTreesAppendString(String filename,String appendString )
    {
    	return createAppendString(filename, appendString ,"trees");
    }
    
    public static void saveLogs(String appendString, ByteArrayOutputStream ess, ByteArrayOutputStream cess, ByteArrayOutputStream weightsStream) throws IOException
    {	    
    	// save ESS, CESS, and normalised weights array
    	{
	    	  OutputStream outputStream = new FileOutputStream(createTxtAppendString("ESS",appendString));
		    	ess.writeTo(outputStream);
    	}

    	{
	    	OutputStream outputStream = new FileOutputStream(createTxtAppendString("CESS",appendString));
			    	cess.writeTo(outputStream);
    	}

    	{
		    OutputStream outputStream = new FileOutputStream(createTxtAppendString("NormalisedWeights",appendString));
		    	weightsStream.writeTo(outputStream);
    	}
    }
    
    private static void saveAvgRej(String appendString, List<Double> avgRejectionList, String paramName) throws IOException {
    	// TODO Auto-generated method stub
	    ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    PrintStream out = new PrintStream(baos);
	    
	    for(int i=0; i<avgRejectionList.size();i++)
		{
	    	out.println(avgRejectionList.get(i));
		}
	    
	    out.close();
	    
	    OutputStream outputStream = new FileOutputStream(createTxtAppendString(paramName,appendString));
	        baos.writeTo(outputStream);
    }

    public static void saveParameter(String appendString, double[] avgRejection, String paramName) throws IOException
    {
	    ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    PrintStream out = new PrintStream(baos);
	    
	    for(int i=0; i<avgRejection.length;i++)
		{
	    	out.println(avgRejection[i]);
		}
	    
	    out.close();
	    
	    OutputStream outputStream = new FileOutputStream(createTxtAppendString(paramName,appendString));
	        baos.writeTo(outputStream);
    }
    
    public static void saveStateSpaceParticles(Sequential[] beastMClist, String appendString, int N_int) throws IOException
    {
		MCMC mccdummy=beastMClist[0].m_mcmc;
		int statenodesNr=mccdummy.getState().getNrOfStateNodes();

    	for(int stateNodeCnt=0; stateNodeCnt<statenodesNr;stateNodeCnt++)
    	{
	    	ByteArrayOutputStream baos = new ByteArrayOutputStream();
		    PrintStream out = new PrintStream(baos);
		    
		    // init the logger with a header
			MCMC mcc=beastMClist[0].m_mcmc;
			String id=mcc.getState().stateNode[stateNodeCnt].getID();
			int index=id.indexOf('.');
			id=id.substring(0,index!=-1?index:id.length()-1);
			mcc.getState().stateNode[stateNodeCnt].init(out);            
		    out.println();
			// log all the particles
		    for(int i=0; i<N_int;i++)
			{
		     	MCMC mc=beastMClist[i].m_mcmc;
		    	mc.getState().stateNode[stateNodeCnt].log(i, out);
		    	out.println();
		    	//System.out.println(mc.getState().stateNode[treepositionInStateArray].toString());
			}
		    
		    out.close();
		    String appendStr;
		    if(id.startsWith("Tree"))
		    {
		    	appendStr=createTreesAppendString(id,appendString);
		    }
		    else
		    {
		    	appendStr=createTxtAppendString(id,appendString);
		    }
		    
		    OutputStream outputStream = new FileOutputStream(appendStr);
		        baos.writeTo(outputStream);
    	}
    }
    
    /*
    public static void saveTreeParticles(Sequential[] beastMClist, String appendString, int treepositionInStateArray, int N_int) throws IOException
    {
	    ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    PrintStream out = new PrintStream(baos);
	    
	    // init the logger with a header
		MCMC mcc=(MCMC)beastMClist[0].m_runnable;
		mcc.getState().stateNode[treepositionInStateArray].init(out);            
	    out.println();
		// log all the particles
	    for(int i=0; i<N_int;i++)
		{
	     	MCMC mc=(MCMC)beastMClist[i].m_runnable;
	    	mc.getState().stateNode[treepositionInStateArray].log(i, out);
	    	out.println();
	    	//System.out.println(mc.getState().stateNode[treepositionInStateArray].toString());
		}
	    
	    out.close();
	    
	    OutputStream outputStream = new FileOutputStream(createTreesAppendString("Tree",appendString));
	        baos.writeTo(outputStream);
	        
    }
    */
    
    static int getCESSexponent_revised( Sequential[] beastMClist, double[] logIncrementalWeights, double[] logWeightsNormalized, int previousExponent, int maxExponent, double stepSize, double normalisedTargetCESS)
    {
    	int nextExponent=maxExponent;
    	double acceptableRange=0.4;
     	int maxcNt=40;
     	double CESSval;
	    double range;
	    double previousExponentDouble=previousExponent*stepSize;
		range=(maxExponent-previousExponent);
     	// here calculate nextExponent
		boolean success=false;
		int oldNextExponent;
		int upperVal=maxExponent, lowerVal=previousExponent;
		//if(previousExponent==0)
		//	previousExponent=1;
		do
		{
			oldNextExponent=nextExponent;
	     	// calculate weights
	     	calculateIncrementalWeights(beastMClist, logIncrementalWeights, previousExponentDouble, (double)(nextExponent*stepSize));
	     	// calculate CESS
			CESSval=CESS_Normalised(logIncrementalWeights, logWeightsNormalized);
	     	// if CESS val == normalisedTargetCESS

			if((CESSval >= (normalisedTargetCESS-acceptableRange)) && (CESSval <= (normalisedTargetCESS+acceptableRange)))
			{
				success=true;
			}
			else
			{
				range=range/2.0;
				if((range<1) && (nextExponent>(previousExponent+1)))
					range=1;
					
				if(CESSval > normalisedTargetCESS)
				{
					// bisection: if CESS val is greater than target, go further
					//upperVal=
					nextExponent=nextExponent+(int)range;
				}
				else
				{
					// if CESSval less than target, then get closer if possible
					nextExponent=nextExponent-(int)range;
				}
			}
						
			maxcNt--;
		}while((success!=true) && (maxcNt>0) && (oldNextExponent!=nextExponent) && (nextExponent>previousExponent) && (nextExponent<maxExponent));
		
		// stop either reached 1 or we came back to previousExponent and still not found or maxCnt expired (should never happen)

		
		if(!success)
		{
			return previousExponent+1;
		}
		else		
    	return nextExponent;
    }

    static double getCESSexponent_double( Sequential[] beastMClist, double[] logIncrementalWeights, double[] logWeightsNormalized, double prevExponent, double normalisedTargetCESS)
    {
    	final double acceptableRange=0.0;//0.1;
     	int maxcNt=40;
     	double CESSval;
	    double range=(1.0-prevExponent);
	    double nextExp=1.0;
     	// here calculate nextExponent
		boolean success=false;
		do
		{
	     	// calculate weights
	     	calculateIncrementalWeights(beastMClist, logIncrementalWeights, prevExponent, nextExp);
	     	// calculate CESS
			CESSval=CESS_Normalised(logIncrementalWeights, logWeightsNormalized);
	     	// if CESS val == normalisedTargetCESS

			if((CESSval >= (normalisedTargetCESS-acceptableRange)) && (CESSval <= (normalisedTargetCESS+acceptableRange)))
			{
				success=true;
			}
			else
			{
				range=range/2.0;
					
				if(CESSval > normalisedTargetCESS)
				{
					// bisection: if CESS val is greater than target, go further
					//upperVal=
					nextExp=nextExp+range;
				}
				else
				{
					// if CESSval less than target, then get closer if possible
					nextExp=nextExp-range;
				}
			}
						
			maxcNt--;
		}while((success!=true) && (maxcNt>0));
		
		
		return nextExp;
    }
    
    static int getCESSexponent( Sequential[] beastMClist, double[] logIncrementalWeights, double[] logWeightsNormalized, double stepSize, double previousExponent, int exponentCnt, int maxvalcnt, double outNextExponent)
    {
   	  int next_exponentCnt=exponentCnt+1;
   	  int upperRange=maxvalcnt;
   	  int lowerRange=next_exponentCnt;
   	  double rangeSize=(upperRange-lowerRange)/2;
   	  double oldCESSval=0.0, CESSval;
   	  double nextExponent;
   	  int direction = 1;
   	  double halfRange=0.001;

   	  
   	  int maxcNt=100;
   	  do
   	  {
   		   nextExponent=next_exponentCnt*stepSize;
   		    //reweight done below, calculation of the incremental part
   	       calculateIncrementalWeights(beastMClist, logIncrementalWeights, previousExponent, nextExponent);
   		   CESSval=CESS_Normalised(logIncrementalWeights, logWeightsNormalized);
   		   
   		   double difference=CESSval-0.9;
   		   
   		   if(difference >=-halfRange && difference <=halfRange)
   		   {// we found it
   		      outNextExponent=nextExponent;
   			   break;
   		   }
   		   
   	       //lowerRange=next_exponentCnt;
   		   if(difference>0)
   		   { 
   			   direction=1;
   		   }
   		   else
   		   {
   		     direction=-1;
   		   }
   	       rangeSize/=2; //=(upperRange-next_exponentCnt)/2;
   		   next_exponentCnt=next_exponentCnt+(int)(direction*rangeSize);
   		   
   		   oldCESSval=CESSval;
   		   maxcNt--;
   	  } while(maxcNt>0 && next_exponentCnt>1); // avoid infinite loops
   	  
	 exponentCnt=next_exponentCnt-1;
   	 return exponentCnt;
   }
    
/* static int getCESSexponent( Sequential[] beastMClist, double[] logIncrementalWeights, double[] logWeightsNormalized, double stepSize, double previousExponent, int exponentCnt, int maxvalcnt, double outNextExponent)
 {
	  int next_exponentCnt=exponentCnt+1;
	  int upperRange=maxvalcnt;
	  int lowerRange=next_exponentCnt;
	  double rangeSize=(upperRange-lowerRange)/2;
	  double oldCESSval=0.0, CESSval;
	  double nextExponent;
	  int direction = 1;
	  double halfRange=0.001;

	  
	  int maxcNt=20;
	  do
	  {
		   nextExponent=next_exponentCnt*stepSize;
		    //reweight done below, calculation of the incremental part
	       calculateIncrementalWeights(beastMClist, logIncrementalWeights, previousExponent, nextExponent);
		   CESSval=CESS_Normalised(logIncrementalWeights, logWeightsNormalized);
		   
		   double difference=CESSval-0.9;
		   
		   if(difference >=-halfRange && difference <=halfRange)
		   {// we found it
		      exponentCnt=next_exponentCnt-1;
		      outNextExponent=nextExponent;
			   break;
		   }
		   
	       lowerRange=next_exponentCnt;
	       rangeSize/=2;
		   if(difference>0)
		   { 
			   direction=1;
		   }
		   else
		   {
		     direction=-1;
		   }
		   next_exponentCnt=lowerRange+(int)(direction*rangeSize);
		   
		   oldCESSval=CESSval;
		   maxcNt--;
	  } while(maxcNt>0 && ); // avoid infinite loops
	  
	 return exponentCnt;
}
*/	
    /*
     * the following is the density of the distribution for the height selection, after formula (24) of the
     * paper on transformations on SMC, log value returned
     */
    public static double calcLogHeightSelection(double theta, double my_height, double meanOfGaussian, double sdOfGaussian)
    {
    	double result =0.0;
    	double fourthird=(4.0/3.0)*theta*my_height;
    	
    	double result1=Math.log(2*Math.abs(theta))-fourthird;
    	double expval=(1-Math.exp(-fourthird));
    	double result2=-(Math.log(Math.sqrt(3)) + Math.log(Math.sqrt(expval)) + Math.log(Math.sqrt(1.0-(3.0*0.25*expval))));
    	double result3=-(Math.log(Math.sqrt(Math.PI*2))+Math.log(sdOfGaussian));
    	double result4temp=2*Math.asin(Math.sqrt(3.0*0.25*expval))-meanOfGaussian;
    	double result4=-result4temp*result4temp/(2*sdOfGaussian*sdOfGaussian);
    	
    	return result1+result2+result3+result4;
    }

    /*    
     * the function addSequence adds a sequence of DNA to the exixting taxon set
     * the return value is the array of the distances number of taxa
     * the new leaf position in the array of leaves can be inferred by the length of the array
     * i.e. distances.length=4 for example, means that the new leaf is t5, i.e. the 5th taxon
    */    

    private static boolean checkAndAddSequenceAndCalculateDistances(Sequential[] beastMClist, int treepositionInStateArray, final String DNAsequence, int nrOfSequencessBeforeUpdate, final int [] distances)
	{
		// add the new sequence first to first element of the list, then use the same shared input for all
		Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
		// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
		// int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
		int numberOfSequencesAfterUpdate=nrOfSequencessBeforeUpdate+1;
		final String seqstr=DNAsequence;
		final TaxonSet txs=tret.m_taxonset.get();
		final Input<TaxonSet> txset=tret.m_taxonset;
		Alignment ali=txs.alignmentInput.get();

		// final int [] distances=new int[nrOfSequencessBeforeUpdate];
		// get the sequence of the last leaf
		final int nrOfElements=seqstr.length();
		
		Arrays.parallelSetAll(distances, e ->
       	{ 
       		String curSeqString=ali.getSequence(e).dataInput.get();
       		int dist=0;
       		// calculate the distance
       		for(int i=0;i<nrOfElements;i++)
       		{
       			if(seqstr.charAt(i)!=curSeqString.charAt(i))
       			{
       				dist++;
       			}
       		}
       		return dist;
       	});
		
		// check if there is any distance equal to 0 (implies that this sequence was already in)
		boolean sequenceAlreadyPresent = Arrays.stream(distances).parallel().anyMatch(i -> i==0);
		
		if(!sequenceAlreadyPresent)
		{// add the sequence to the alignment
			String taxonID="t"+numberOfSequencesAfterUpdate;
			Sequence seq=new Sequence(taxonID, seqstr);
			
			ali.initializeWithAddedSequenceList(Arrays.asList(seq), false);
		}
		
		// return true if the sequence was added
		return (!sequenceAlreadyPresent);
	}

    /*    
     * the function addSequence adds a sequence of DNA to the exixting taxon set
     * the return value is the array of the distances number of taxa
     * the new leaf position in the array of leaves can be inferred by the length of the array
     * i.e. distances.length=4 for example, means that the new leaf is t5, i.e. the 5th taxon
    */    

    private static ArrayList<List<Double>> drawLeafAndHeightAndCalculateLog(Sequential[] beastMClist, int treepositionInStateArray, int popsizepositionInStateArray, int sequenceLength, int nrOfSequencessBeforeUpdate, final int [] distances)
	{
		// add the new sequence first to first element of the list, then use the same shared input for all
		Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
		// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
		// int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
		int numberOfSequencesAfterUpdate=nrOfSequencessBeforeUpdate+1;
		
		final TaxonSet txs=tret.m_taxonset.get();
		final Input<TaxonSet> txset=tret.m_taxonset;
		Alignment ali=txs.alignmentInput.get();


		final Input<Alignment> aliinput=txs.alignmentInput;
		

		final Integer lengthOfSeq=sequenceLength;
		final double lenSeq=lengthOfSeq.doubleValue();
		
		// needed for drawing of the height
		final double sdOfGaussian=1.0/Math.sqrt(lenSeq);

		
		// after having updated the tree we need to update all obj that hv the tree as input?
		// probably not needed as the initialisation is done anyway in the mcmc init
		
		// the first element is the leaf, the second is the distance of the new leaf from each of the previous ones
		ArrayList<List<Double>>  retValues = new ArrayList<>();
		for(int i=0; i< beastMClist.length; i++)
		{// initialise and add a list of doubles for
			retValues.add(new ArrayList<Double>());
		}

		Arrays.parallelSetAll(beastMClist, e ->
	       	{ 
	       		// add sequence here to all
	       		// would be good to have a shared taxa, taxonset, alignment for all particles
	       		
	       		List<Double> retvalRow=new ArrayList<>(Arrays.asList(new Double[] {0.0,0.0,0.0,0.0}));
	       		MCMC mc=beastMClist[e].m_mcmc;
	       		State stt=mc.getState();
	       		double theta=stt.stateNode[popsizepositionInStateArray].getArrayValue();

	       		/* start of the part of drawing the leaf */
				Integer selectedLeaf=0;
				//double theta=beastMClist[e];
				{
					double qtToElevate=(lenSeq*theta)/(nrOfSequencessBeforeUpdate+(lenSeq*theta)); // this will hv to change and contain the effective pop size
					double []probabilityWeight= new double[distances.length];//{0.1,0.2,0.25,0.85};//
					Arrays.parallelSetAll(probabilityWeight, ee -> {return Math.exp(distances[ee]*Math.log(qtToElevate));});
					org.apache.commons.math3.util.Pair<Integer, Double> itemToInit=new org.apache.commons.math3.util.Pair<Integer, Double>(0,0.0);
					List<org.apache.commons.math3.util.Pair<Integer, Double>> pmfWeights=new ArrayList<org.apache.commons.math3.util.Pair<Integer, Double>>(Collections.nCopies(probabilityWeight.length, itemToInit));
	        		
					Arrays.parallelSetAll(probabilityWeight, ee ->{
	        			pmfWeights.set(ee, new org.apache.commons.math3.util.Pair<Integer, Double>(ee,probabilityWeight[ee]));
	        			return probabilityWeight[ee];
	        		});
					EnumeratedDistribution enDist=new EnumeratedDistribution<>(pmfWeights);
					// the following is the draw of the leaf (position of the leaf in the array)
					selectedLeaf=(Integer) enDist.sample();
					//selectedLeafMap.replace(e, selectedLeaf);
				}
	       		/* end of draw of the leaf */
	       		double my_height=0.0;
				double meanOfGaussian=2*Math.asin(Math.sqrt(distances[selectedLeaf.intValue()]/lenSeq));
	       		/* start of the part of drawing the height */
				{
					double upperBoundTruncatedGaussian=2.094395102393195; // this is 2*Math.asin(sqrt(3)/2.0);
					//org.apache.commons.math3.distribution.NormalDistribution norm=new org.apache.commons.math3.distribution.NormalDistribution(meanOfGaussian, sdOfGaussian);
					//TruncatedNormal tn = new TruncatedNormal(meanOfGaussian, sdOfGaussian,  Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian);

					double inerval;
						double beta=TruncatedNormal.sampleUpgraded(Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian); //tn.sample();
						if(beta > upperBoundTruncatedGaussian)
						{
				            throw new RuntimeException(
				                    "Unable to draw properly from the truncated Gaussian\n");
						}
						// from the paper on Sequential Monte Carlo transformations,
						// calculate the height (formula 24 of paper)
						// decide if it is better to have a different height for every particle:
						// the process of selectin height is independent of the tree
						double sinVal=Math.sin(beta/2.0);
						inerval=1.0-((4.0/3.0)*sinVal*sinVal);

						if(inerval > 0)
						{
						    my_height=(-3.0/((4.0)*theta))*Math.log(1.0-((4.0/3.0)*sinVal*sinVal));
						}
						else
						{
				            throw new RuntimeException(
				                    "Error in formula from the truncated Gaussian, negative log argument!!!\n");
						}
				}

				/* end of draw of the height */
				retvalRow.set(selectedLeafPosition, selectedLeaf.doubleValue());
				retvalRow.set(heightPosition, my_height);
				retvalRow.set(meanGaussianPosition,meanOfGaussian);
				retvalRow.set(sigmaGaussianPosition,meanOfGaussian);

				/* 
				 * here we can calculate the components for the weights
				 */
				/*
				 * Leaf first
				 */
	       		   double Ms=distances[selectedLeaf.intValue()];
	       		   double N=lenSeq; // length of sequence
	       		   double t=nrOfSequencessBeforeUpdate;
				   double logUnnormalisedIncrementalWeightsLeafSelection=Ms*(Math.log(N*theta)-Math.log(t+N*theta));

				  // retvalRow.set(logLeafComponentPosition, logUnnormalisedIncrementalWeightsLeafSelection);
				   /*
				 * height afterwards
				 */
				   double logUnnormalisedIncrementalWeightsHeightSelection=calcLogHeightSelection(theta, my_height, meanOfGaussian, sdOfGaussian);

				  // retvalRow.set(logHeightComponentPosition, logUnnormalisedIncrementalWeightsHeightSelection);
				   
				   retValues.set(e, retvalRow);
				  return beastMClist[e];
	       	});
		
	   return retValues;    	
	}

    private static double calculateLeavesComponents(Sequential[] beastMClist, int treepositionInStateArray, int popsizepositionInStateArray, int sequenceLength, int nrOfSequencessBeforeUpdate, final int [] distances, ArrayList<HashSet<Integer>> leaves)
	{
    	double result =0.0;
		/* 
		 * here we can calculate the components for the weights
		 */
		/*
		 * Leaf first
		 */
    	IntStream.range(0, 10).parallel().forEach(i -> {
    		// outer loop if for each particle
    		MCMC mc=beastMClist[i].m_mcmc;
       		State stt=mc.getState();
       		double theta=stt.stateNode[popsizepositionInStateArray].getArrayValue();
    		double N=sequenceLength; // length of sequence
       		double t=nrOfSequencessBeforeUpdate;
       		HashSet<Integer> leavesSet=leaves.get(i);
       		for(Integer selectedLeaf:leavesSet)
       		{
    			double meanOfGaussian=2*Math.asin(Math.sqrt(distances[selectedLeaf.intValue()]/N));
       			double Ms=distances[selectedLeaf.intValue()];
       			double logUnnormalisedIncrementalWeightsLeafSelection=Ms*(Math.log(N*theta)-Math.log(t+N*theta));
       			/*
       			 * height afterwards
       			 */
    		   // double logUnnormalisedIncrementalWeightsHeightSelection=calcLogHeightSelection(theta, my_height, meanOfGaussian, sdOfGaussian);
       		}
    	  });

    	return result;
	}
    
    private static ArrayList<List<Double>> drawLeafAndHeightAndCalculateLogUpdated(Sequential[] beastMClist, int treepositionInStateArray, int popsizepositionInStateArray, int sequenceLength, int nrOfSequencessBeforeUpdate, final int [] distances, ArrayList<HashSet<Integer>> leaves)
	{
		// add the new sequence first to first element of the list, then use the same shared input for all
		Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
		// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
		// int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
		int numberOfSequencesAfterUpdate=nrOfSequencessBeforeUpdate+1;
		
		final TaxonSet txs=tret.m_taxonset.get();
		final Input<TaxonSet> txset=tret.m_taxonset;
		Alignment ali=txs.alignmentInput.get();


		final Input<Alignment> aliinput=txs.alignmentInput;
		

		final Integer lengthOfSeq=sequenceLength;
		final double lenSeq=lengthOfSeq.doubleValue();
		
		// needed for drawing of the height
		final double sdOfGaussian=1.0/Math.sqrt(lenSeq);

		
		// after having updated the tree we need to update all obj that hv the tree as input?
		// probably not needed as the initialisation is done anyway in the mcmc init
		
		// the first element is the leaf, the second is the distance of the new leaf from each of the previous ones
		ArrayList<List<Double>>  retValues = new ArrayList<>();
		for(int i=0; i< beastMClist.length; i++)
		{// initialise and add a list of doubles for
			retValues.add(new ArrayList<Double>());
		}

		Arrays.parallelSetAll(beastMClist, e ->
	       	{ 
	       		// add sequence here to all
	       		// would be good to have a shared taxa, taxonset, alignment for all particles
	       		
	       		List<Double> retvalRow=new ArrayList<>(Arrays.asList(new Double[] {0.0,0.0,0.0,0.0}));
	       		MCMC mc=beastMClist[e].m_mcmc;
	       		State stt=mc.getState();
	       		double theta=stt.stateNode[popsizepositionInStateArray].getArrayValue();

	       		/* start of the part of drawing the leaf */
				Integer selectedLeaf=0;
				//double theta=beastMClist[e];
				{
					double qtToElevate=(lenSeq*theta)/(nrOfSequencessBeforeUpdate+(lenSeq*theta)); // this will hv to change and contain the effective pop size
					double []probabilityWeight= new double[distances.length];//{0.1,0.2,0.25,0.85};//
					Arrays.parallelSetAll(probabilityWeight, ee -> {return Math.exp(distances[ee]*Math.log(qtToElevate));});
					org.apache.commons.math3.util.Pair<Integer, Double> itemToInit=new org.apache.commons.math3.util.Pair<Integer, Double>(0,0.0);
					List<org.apache.commons.math3.util.Pair<Integer, Double>> pmfWeights=new ArrayList<org.apache.commons.math3.util.Pair<Integer, Double>>(Collections.nCopies(probabilityWeight.length, itemToInit));
	        		
					Arrays.parallelSetAll(probabilityWeight, ee ->{
	        			pmfWeights.set(ee, new org.apache.commons.math3.util.Pair<Integer, Double>(ee,probabilityWeight[ee]));
	        			return probabilityWeight[ee];
	        		});
					EnumeratedDistribution enDist=new EnumeratedDistribution<>(pmfWeights);
					// the following is the draw of the leaf (position of the leaf in the array)
					selectedLeaf=(Integer) enDist.sample();
					//selectedLeafMap.replace(e, selectedLeaf);
				}
	       		/* end of draw of the leaf */
	       		double my_height=0.0;
				double meanOfGaussian=2*Math.asin(Math.sqrt(distances[selectedLeaf.intValue()]/lenSeq));
	       		/* start of the part of drawing the height */
				{
					double upperBoundTruncatedGaussian=2.094395102393195; // this is 2*Math.asin(sqrt(3)/2.0);
					//org.apache.commons.math3.distribution.NormalDistribution norm=new org.apache.commons.math3.distribution.NormalDistribution(meanOfGaussian, sdOfGaussian);
					//TruncatedNormal tn = new TruncatedNormal(meanOfGaussian, sdOfGaussian,  Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian);

					double inerval;
						double beta=TruncatedNormal.sampleUpgraded(Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian); //tn.sample();
						if(beta > upperBoundTruncatedGaussian)
						{
				            throw new RuntimeException(
				                    "Unable to draw properly from the truncated Gaussian\n");
						}
						// from the paper on Sequential Monte Carlo transformations,
						// calculate the height (formula 24 of paper)
						// decide if it is better to have a different height for every particle:
						// the process of selectin height is independent of the tree
						double sinVal=Math.sin(beta/2.0);
						inerval=1.0-((4.0/3.0)*sinVal*sinVal);

						if(inerval > 0)
						{
						    my_height=(-3.0/((4.0)*theta))*Math.log(1.0-((4.0/3.0)*sinVal*sinVal));
						}
						else
						{
				            throw new RuntimeException(
				                    "Error in formula from the truncated Gaussian, negative log argument!!!\n");
						}
				}

				/* end of draw of the height */
				retvalRow.set(selectedLeafPosition, selectedLeaf.doubleValue());
				retvalRow.set(heightPosition, my_height);
				
				// we need to save mean of gaussian and stddev of gaussian
				Tree particleTree = (Tree)stt.stateNode[treepositionInStateArray];
				int chosenNode=RandomTree.getClosestFromLeaf(particleTree, selectedLeaf.intValue(), my_height);
				HashSet<Integer> leavesSet=leaves.get(e);
				leavesSet.add(selectedLeaf);
				leavesSet=RandomTree.FindAllLeavesFromInternalNode(particleTree, chosenNode, leavesSet);
				leaves.set(e, leavesSet);
				//retvalRow.set(index, (double)chosenNode);
				/* 
				 * here we can calculate the components for the weights
				 */
				/*
				 * Leaf first
				 */
	       		   double Ms=distances[selectedLeaf.intValue()];
	       		   double N=lenSeq; // length of sequence
	       		   double t=nrOfSequencessBeforeUpdate;
				   double logUnnormalisedIncrementalWeightsLeafSelection=Ms*(Math.log(N*theta)-Math.log(t+N*theta));

				  // retvalRow.set(logLeafComponentPosition, logUnnormalisedIncrementalWeightsLeafSelection);
				   /*
				 * height afterwards
				 */
				   double logUnnormalisedIncrementalWeightsHeightSelection=calcLogHeightSelection(theta, my_height, meanOfGaussian, sdOfGaussian);

				  // retvalRow.set(logHeightComponentPosition, logUnnormalisedIncrementalWeightsHeightSelection);
				   
				   retValues.set(e, retvalRow);
				  return beastMClist[e];
	       	});
		
	   return retValues;    	
	}

    /*
     * this function is used to calculate the adaptive version of the SMC
     */
    public static double calculateWeightedMoment(Sequential[] beastMClist, double[]logNormalisedWeights, int positionInStateArray, int momentNumber)
    {
    
    	if(momentNumber<1 || momentNumber>2)
    	{
            throw new RuntimeException(
                    "The function does not calculate moments higher than 2\n");
    	}
    	
    	double momentSum=0.0;
    	double weight;
    	double paramVal;
    	MCMC mc;
    	int nParticles=beastMClist.length;
    	for(int i=0; i<nParticles; i++)
    	{
    		weight=Math.exp(logNormalisedWeights[i]);
    		paramVal=beastMClist[i].m_mcmc.getState().stateNode[positionInStateArray].getArrayValue();
    		
    		if(momentNumber==1)
    		{
    			momentSum+=(weight*paramVal);
    		}
    		else
    		{ // assume it's second moment
    			momentSum+=(weight*weight*paramVal*paramVal);
    		}
    	}
    	
    	momentSum/=(nParticles*1.0);
    	return momentSum;
    }
    
    public static double calculateParameterVariance(Sequential[] beastMClist, double[]logNormalisedWeights, int positionInStateArray)
    {
    	double firstMoment=calculateWeightedMoment(beastMClist, logNormalisedWeights, positionInStateArray,1);
    	double secondMoment=calculateWeightedMoment(beastMClist, logNormalisedWeights, positionInStateArray,2);
    	
    	return secondMoment-(firstMoment*firstMoment);
    }

    /*    
     * the function addSequence adds a sequence of DNA to the exixting taxon set
     * the return value is the array of the distances number of taxa
     * the new leaf position in the array of leaves can be inferred by the length of the array
     * i.e. distances.length=4 for example, means that the new leaf is t5, i.e. the 5th taxon
    */    

    private static HashMap<Integer, Integer> addSequence(Sequential[] beastMClist, int treepositionInStateArray, int popsizepositionInStateArray, int seqNr, int nrOfSequencessBeforeUpdate, final int [] distances)
	{
		// add the new sequence first to first element of the list, then use the same shared input for all
		Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
		// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
		// int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
		int numberOfSequencesAfterUpdate=nrOfSequencessBeforeUpdate+1;
		String taxonID="t"+numberOfSequencesAfterUpdate;

		final String[] sequencesArray = {
				"GTTGGCACAGTCGAATGACTGGTATACTGTTCGTCAACGATTACATAGGACTCGACTGAGCGGGACGAACTTAGGCATAATGGGGAAAGTAGCCTCCCTTCCATACCGCAAGATTTGGTATACTCTCCCCGTCTGGCAGAATCGTCCCCCTATTTTGTCCACTAGTTTACACGTGGCAGAAGCCGGACGGGGTATTGCCTCGTCTCGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAAACAAGGACCACTACAGATGTGAATGGGCACCACTTAAGCATCTGCATCGAACACATGGGAAAACCTTTCGTTCAAGCAAATCTAAACTTAGACCACCGACCTCTTGTGATGCTCATTCGGACGGAGAGTGATCCAAGGGAGTCTGGGTTACGGGTCTTGTGTAGACCTTTCAGTCAATGCCTCCATTTTGATCCAAACACTGGATGTTCAATCACTCGTTAGTAGCCCACTGTTATGAGAGGACACCAGTTGACCAGGATGCCGACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGTTGGCAAACTATCGAGACTGTTTGCACGTAAAACCGGCATTCAGCAACTGTTGACCGGTGGTAATCCTCGAGGGGTCAAATCACTGATCCATTAGAGTACCCTGTACTAGCCTATACAACCGAAGGTAAATCGATACCTAACAGGGGTTATGCGCTCCTAAACGCTTCCCAGAGTGTGCGCTGCTCGGCTAAGGCGTTCCAAACTTGTAAAAATCTTTACGCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTGTACGTTGATCGGAGCTGCGATTAACACACGGTAAGGCACGGTGGTATCAAAACTTTACATTCATGAAGTATGAGGGTGTCACAATCAAAGTACCAGGCATAAAAATTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACAAATTCGGGCGGAAGACCCTAATTTACAAAACCTCGGAGGTAAG",
                "CTGGGCTCTGGGGGTATTCTACTCGGCGCTCCGTTCGCACATCAACTGATAATAGTATCATCAGCTTACGGGTCGGCGTGGGCGGTCCATCTTGTTTGCCTAAAGTTGATAAAAGGAGGTGGCAGATGCATATTGTCGAAGGAATCAGCATGGTCTGGCATTGAAATCCAAGCCTTTATGGGCCTGTGCTGCTCCGGACCACTTGCCTCTGCGCGGATCCAGCTCGTCACTGAGACTCTCCACCAATAATTCAAGCAGCTTTGAACGTGGTAAGGATAGGACCCCGTCCCGTGGCTACGCAGTCTACTCGGACTCTATGAAGCGAAAATCAGGTTCTCATGAAGGTCTCAAATCGATGTTTCTCAATGACTATGCAATAGATGTGCGGTGTACTACGCTTCAGATGATTGAATACTTGCCCTGCATCGAGACGCAATGTGTGGTGATATGATAAAAATCAGGTCCATGGAACTTCCAGAATAGTCAACGACTGTATGCGGCGCGTCAATGAAAAACCGACTTCCGGACTAACGTGCGATTCGCCAAGGACCACACCCGAGATGGCGCACACAGAATAGACTGGGCAAATATATTATGCTACTTTTGGTCATCGGGGTACGAGAGGTAGCCTCAGAACCGGATAAGCAGCCGCCCTGCCGTAGGGGGTCTCCGTGACAATAGACTGTAATCATCACAGTCGTAATCAGGCGTTCCATACAGTTATGCTTCGCTGAGGGTCTGGTAGATGGCTTCGGACTAGACGGCAGCACCGTGTTGACGGCCTCATTACGGGTGCAAGACCGGTTTGAGCTTCACGTGTCACAGATTTTTAAGTTGCAAATCACTCATCTCCGACACAGAGGGAAGAGGTAAGCGCACTGTTCCTCTCTGACTAGACTGGAAAGGGGTTAAAGCACTTTCTCTATTGGCTCCCATATCCGTCACTTCTGACTTCATCATGATCTACACCGCAATGCCCACTTTTCTGAAAGGACCTTGG",
                "CTGGGCTAGGGGGCTATTCTACTCGGCGTTCCGTTCGCACATCAACTGATAATTGTTTCATCAGCTTACGGGTTGGTGTGGGCGGTCCATCTGGTTTGCCTAAAGGTGATAAAAGAAGTTGCCCGATGCATATTGTTGAAGCTATCATCATGGTCTGGCATTGCAATCCAAGCCTTTATGGGCCTGTGCTGCTCCGGCCCACTTGCCTCTGCGCGGATCCAGCTCGTGAATGAGACTCTCCTCCAAGAAATCAAGCAGCTCTGAACGTGGTAAGAATTGGAACACGTCCCGTAGCTACGCAGTCTACGCGGACTCTATGAAGCGAAAATCAGGTTCCCATGAAGGTCTCTAATCCGTGTTTCTCAATGAATTTGCAATAGATGTCCGGTGTACTACGCTTAAGATGATTGAAAACCTGCCCTGATTCGAGATGAAATCTGTGGTGACAAGATAGACGTTGAGTCCCTGGAACTCCCAGAATAGTCAACGACTGTATGCGGCGCGCCGATCAAAACCCGACATCCGGATTGAGGTGCGATTCGCCAAGGACCCCACCCGAGATGGGGCACACACAATATACTGGGTTAATATATTGTGCTACTTTTGGTCAACGGGGTACCAGAGGGAGCCTCTGAACCGGATAAGCAGCCCCCCCGCCGTAGGGGGTCTCCGTGACAATAGACTGTAATCATCACATTTATATTCAGGCGTACCATAGAGTTATGCTTCGCTGAGGGTCCGGTAGATGGCTTCTGTCTAGACGGCAGCAGCGTGTTGACGGCCTCATAACGGTTGCAAGACCGGTTTGAGCTCCACGGGTCACAGATTTTTAAGTTGCAAATTGCTCATCTCCGACAGAGAGGGAAGAGCTAAGGGCGCTGTCACTCTCTGACCAGACTGGAAAGGGTTTAAAGCACTGTCTCTATTGGCTCCAATCTCCGTCACTTCTGACTTCATCATGATCTACACCGCAATGCCCAGTTTTCAGAAAGGAACTTCG",
                "AAGCCCGTGCTTTACAGTCCGCATTTTACTAGTGCACTAAATACTACACTCCGTTGGAGTTCCGCCGGATAGGTATGCATCTAAAAGAGTCAGGCCTCCATTCTCTCAGCAACGCGGATTCCAGCAGAGTCTCCATTTCCGATTGTCGTCCTGATACCGGGGGAAGAGGTGTTATTTATCTCTCCACCATATTCTTGGTTGTACCCTGTCTTGGCCTTAGGGGAGCCTTACTTGCTCCTGTGACCGGGTGAACTCCCGCGTTCCCGTCCGAGGGTCTTCCTTCGAAACACGCTGTTCGTCACCACGTTCGTCCGAGAGACCCGTGCATTTCCAGACGTACTCGTGTGCTCCACGTCAATCACGGACTAATTGTAGTACAGCGCGTGTTAATTGCAGCCTTCTAAGATCCGTTAGCCAGGGGATGGATAAATCCCCTAGCGTTAGTTGAACTCAGTAGAAGACGGAACTACTACCTATATCGCTACCGACACCGGGCAGGGGGACAGATGGGGACTAGAGCCTTATATCGTAGGATGAAAAGTCCTCCCAAGACATTGTCAGACGGATCCCAGTCCATTGATAGGTGTAGCGGGTAGTCATCTAATGTGTAAGGCCAACTATGATAGTACACATCCACGCACATACGCTTCGTAGATGCCGGCCTGCCTCCTCAATCTAGTAAGGATCTACTGCATTTTATATAACAAACAACGACGGACTTGTTCCGTGCTACTATACAGTCTGGAACACAGCCATGCGTGGAAAACTCACCAGCCAGCATGGGTGTAAGGACTTCTATAGGCAGCGGTGGAAGAGTAGGGTAGTTGGTATTTCGTCAGTTGGCAAGGTATGTTAGCGGGGCGTCGAGGGTTATGACGGTTGCACAAGTCTGGGTGATTATCAACAGCAAGCGTCGTTGACCAGTACGTTCGATACCGGGAAGGTCCACGTGCGTTTACACAATGAGACTATAACCCGCGCCAAACGACACAAGAAAATA",
                "AAGACCGTGCTTTGCGGTCCACACTTAAATAGTGCAAGAATTACTACACTCCATTGGTGTTTCGCCGGATATGTATGCATCTAAAAGAGTCAGCCCTCCAGTCTCGCAGCAACCCGGATTCCAGCAGAGTCTCTATTTCGGATTGTCGTCCTGAGACCGAGGGAAGAGCTGCTCTTTATCTCTCCACCATATTCTTGGTTATACCCTGTCTTGTCCTTAGGGGAGCCTTCCTTGCGCCTGTGAACGGGAGAGCTCCCGCGTGCCCGTCCGTGGGTCTACGTTCGAAAAACTCTGTTCGTCACTACGTTCGGCCGAGAGTACCGTGCACTTCCAGACGTACTCGTGTGCTCCACGTCAGTCACTGACTAATTGTAGCACAGCGCCTTTGAATTCAAGCCTTCTAAGATCCGTTAGCCAGGGGATGGACAAATCCCCTAGCGATAGTTGAACTCAGTAGAAGACGGAACTACTACCTATATCGCTACCGACACCGGGCAGGGGGGCAGATGGGGACTACCGCCTTATATCGTAGGATGAAAAGTCCTCCCATTACAATGTCAGAGCGATCCCACCCCATTGATGGATGTAGTGGATAGTCACCTAATGTGTATGGCCAACTATGATAGTACACTTCCACGCACATACGCTTCGTTGATGCCCGCCTCCCACCTAGATCTAGTAAGTATCTCCTGCATTTTATATAACAAACAGTGACGGACTTGTTCCGTGCGACTATACAGCCTGGAAGACAGCCATGCGTGGAAAACTCACCCGCCAGCATGGGTTTAAGGACTTCTATAGGCAGCGGTGGATGAGTAGGGTAGTTGGTATTTCGTCAGTTGGCAAGGTATGTTAGCGGGGCGTCGAGGGTCATGATGGTTGCACAAGTCTCGGTGATTATCAACTCCAACCGTAGTTGACCAGTACGTTCGATACCGGGAAGGTCCAGGTGCGTTTACACAATGAGACTTTAACCCGCGCCAAACGACACAAGAAAATA",
                "AGTGCTCAAGCCGGACCTGACGCGACCAAATATCCATCTTGAGTTCCCAAGTCTCTACACACAGCGGGGAGTTCTCGCATCAACTGACCTATCGTCGCGATTATCTCAGCGGTAACCCCAGCAGTAAGAACCTAGAGATAGTCGCCGTTAAGTTGTACATTATGAGTTATTTGACAAACTTCACAAACTGCAATTCCGGCGGGCCGGACTTTCCCATTGCGCGGCTCTCTACCACGTCTGGGGAAGCACTTACATCAGTAGCTCTTGTGCTTGGCCAACACACATACGATAAAGATCCAACATCTCCGTGCGTGGGGGCAATTCCCACACACAGCATTGCATTGGTTCAGACCAGCATCTCAGAGTGCGAATAAGCGGGCAAATTTTCATTGCTACAACCGCGATCTTCGTTTATGCTCGGCCGGAAATTTGGAAAGGAGCAAAGCTGACCCACGAGCGCGAGTCCCGCTAGCAGAAAAGTCATGTTGCATGCGTAACGGCAGTACGGGCACGGGGGTCGACCGCACTACAGATGTATGCAGTAATATTTGACTAGGGCCCTCAGGTGTGTAAACAGTAAACCGGAAATTCTCTACGTTGTTTTAGTGGACTCCCCTCTCAGGTTAAGGGGCGCCGACGTAACGCGACCGGCTTTAACATTGCGATAATCAATAGGCTGCGCAATTGTAATTCTAGGTTCTAGATAAAAGTTGGATAGTGCACGTTGTAACTACCTGACTATACGCTGCAGCGTCACAAGCATAAGTCCCCTGTGGTAGTGCTCAGTAAGGCTCACTCAGGGTACGTGCAGCGTCTTTTTCGTGCAGCCGAGCATAGTCTAAACGTTTGAGTCTAAACATAGTCAGAACGGTATGCCACTTCCCTCTCGACGACTAGCCACACACCGTGTTACAGGCTGAGTCAAAAGTATTGTGCAGAAACTAAATGGCAGTACCACAAGAGTGCCTTTTTCGGGTTTACTGTGCACTTGCGAGATCGT",
                "AGTGCTCAAGACGGACGTGACGCGACCAAATATCCATCTTGAGTTCCGAAGTCTCTACACATAGCGGGGGGTTCTCGCATCAACTGACCTATCGACGCGACTATCTCTGCGGTAACCCCAGCAGTAAGAAGAAAGAGATAGTCGCCGTTAAGTTGTACATTATGAGTTATTTGACAAACTTCCCAAACTGCAATTCCGGCGGGCCGCACCTTCCCATTGCGCGGCTCTCGACCACGTCGGGGGAAGCACTTACATCAATATCTCACGTGCTTGGCAAACACACATACGATAAAGGTCCAACATCTCGGTGCGTGGGGGCAAGTCCCACACACAGCATTGCATTGGTTCAGACCAGCATCTCAGAGTGCGAATAAGCGGGCAAATTTTCTTTGCTATAACTGCGATCTTCGTGCATGCTCGGCCGGAAATTTGGAAAGGAGCAAGGCTGACCCACCAGCGCGAGTCCGGCTAGTATAAAAGTCATGTTGCATGCGTACCGTCAGTACGGGCACGGGGTTCGACCGCACAACAGGTGTTTGCAGTAATATTTGGCTAGGGCTCTCAGGTGTGTAAACAGTAAGCGGGAAATTCTCTACGTTGTTTAAGAGGACCCCCCTCTCAGGTTGAGGGTCGCCGATGTAACGCGGCCAGCTTTAACATTGCGATATTCAAGTGCCTGCGCAATCGTAATTCTAGGTTCTAGATAAAAGTTGGATAGTGCACGTTGTAACTACCTGTCTATACGCTGCAGCGTCACAAGCATAAGTCCCCAGTGGTAGTGCTCATTAAGGCTCACTCAGGGTCCGTGCAGCGTCTTTTTTGTGCAGCCGAGAATAGTTTAAACGTTTGAGTATAAACATAGTCAGAACGGTATGCCAGTCCCCTCTCGACGAGTAGCCACAGGCCGCGATAAAGGCTGAGTCAAAAGTATTGTCCAGAAACTAAATGACAGTACCACAAGAGTGCCTTTTTCCGGCGTACGGTGCAATAGCGAGATCGT",
                "GCTGGCACAGTCGAATGATTGGTATAATGTTCGTCAACCATTACATAGGACTCGTCTGAGCGGGACGAACTTAGGCATCATGGGGAAAGAAGCCTCCCTTCTATAGCGCAAGATTTGGTATACTCTCCCCGTCGGGCAGAATCGTCCCCCTATTTTTTCCACTAGTTTACAGGTGGCAGAAGCCGGACGGCGTATTGGCCCGACTAGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAGACAAGGACCACTACCGATGTGAATGGGCACCACGTAAGCATCTGAATCGCACACATGGGAAAACCTTTCGTTCAACCAAACCTAAAGTTATACCACGGGCCTTTTGTGATGCTCATTCGGACGGAGGGCCATCCAAGGGAGTCTGGGATGGGGGTCTTGGGTAGACCTTTCAGTCAAGGCCTCCATTTTGATCCAAACATCGAATGTTCAATCACTCGTTAGTAGCCGACTGTTATGAGAGGACACCAGTTGACCCGGATGCCCACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGGTGGCATCGTATCGAGACTGTTTGCACGTAAAATCGGCTTTCAGCGACTGTTGACCGGTGGTAATCCTCGAGGGGTCAATTGACTAATCTATTAGAGTACCCATTACTACTCTATATAACCGAAGGTAAATCGATGCCTAACAGGGGTTATCCGCTCCTAAACGCTTCCCAGAGTGTGCCCAGCTCGACTAAGGCGTACCAAACTTGTTAAAATCTTTACTCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTACACGTTAATCGGAGCTGCGAGTAACACACGGTAAGGCACGGTGGTATCAAACCTTTACATTCAATAAGCATTAGGGTGTCACAATGAGAGTAGCAGGCATAACATCTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACGAATTCGGGCGGAAGACACTAATTTACGACACCTTGGTGGAAAG",
                "AAGACCTTGCTTTGCGGTCCGCACTCTACTAGTGCACTAAATAGTTCACTCCGTGGGAGTGTCCCCGGATATGTATGCATCTAAAAGACTGAAAGCTTCAGTCTCTGAGCAACGCGGCTTCCAGCAGAGTCTCGTTGTCCGATTGCCGTCCTGATACCTAGGGAAGAGTTGCTCTTCCTCTCTCCACCAAATTCTTGGTTGTACCCTGTCTTTGCCTTGGGGGAGCCTTACTTGCCCCTGTGACCGGGTGAACCCACGCGTGCCCGTCCGAGGCTCTACCTACGAAAAACTCTTTTCGTCACCAGCAAGTGCCGAGGGAACCCAGCACTTCCAGACGCACTGGTGTCCACCACGTCAATACCTGACTATTTGGAGTACGGCGCCTATTAGTTTAAGACTTTGAAGATCCGCTAGCCACGGTATGGATCAATCCCCTATCCTCAGTTTAACTCAATAGAAGACGGTTCAACTCCCTATATCGCAACCGACTCCGGGCAGTGGTAGAGATGGGTACTAGCGCCTTATATCGTAGGATAAAAAGTCCTCCCATTACATTATCAGACCGTTCCCAGGCCATTGATAGATGTAGTGGGAATACACCTAATGTGTATGTCCGCCTATGATTGTACAAATCCACGCACATGCTCTGCGTGGATGCCCGCCTCCCTCCCAGATGCAGTAAGGATCTCCTGCATTTTGTTTAACAACCAACTACGGAGTTGTTCCGTGCTACACTTCAGTCTGGAAGCCAGCCATCCATCGAAAAGTCACCAGGCAGCATGGGTGGAAGGCCTTCTATAGGCAGCGGAGGAAAAGTAGGGTAGTTGGTATTTCTTTCGTTGGCAAGGTATGTTAGCGGGACGTCGAGGGCCGGGATGGTTGCACAAGTCTAGGCGATTATCAACTGCAAAAGTTGATGACCAGTCCGTTCGATACCGGGAAGGTCCCTGAGGGGTTACACATCGAGACTAGAACGCGCGCCAAACGACACAAGAAGATA",
                "AAGACCGTGCTTTGCGGTCCGCACTCTACTAGTGCACTAAATAGTTCACTCCGTGGGAGTTTCCCCGGATATGTATGCATCTAAAAGACAGAAACCGTCAGTCTCTCAGCAACGCGTCTTCCAGCAGAGTCTCGTTGTCCGATTGCCGTCCTGATACCTAGGGAAGAGTTGCTCCTCATCTCTCCACCAAATTCTTGGTTGTACCCTGTCTTTGGCTTGGGGGAGCCTTACTTGCTCCTGTGACCGGGTGAACCACCGCGTGCCCGTCCGAGGGTCTACCTTCGCAAAACTCTTTTCGTCACCAGGAACGGCAGAAGGAACCCAGCACTTCCAGACGCACTGGAGTCCTCCACGTCAATCCCTGACTATTTGTAGGACCGCGCCTATTAGTTTAAGACTTTGAAGATCCGCTAGCCACGGTATGTATAATTCCTCTATAGTTAGTTTAACTCAATAGAAGACGGTTCTACTCCCTATATCGCAACCGACTTCGGGCAGGGGAACAGATGGGTACTAGCGCCTTCTATCGTAGGATAAAAAGTCCTCCCATTACATTGTGAGGCCGTTCCCAGGCCATTGATAGATGTAGTGGGAAGACACCTAATGTGTATGATCACCTATGATTGTACAAATCCACGCACATGCTCTGCGTTGATGCCCGCCTCCCTCCCAGATCCAGTAAGGATCTCCTGCATTTTGTTTAACAACCAACTACGGAGTTGTTCCGTGCTACACTTCTGTCTGGAAGACAGCCACCCGTCGAAAAGACAACAGGCAGCATGGGTGGAAGGCCTTCTATAGGTAGCGGAGGAAAAGTAGGGTAGTTGGTATTTCTTTCGTTGGCAAGGTATGTTAGCGGGACGTCGAGGGCCGGCATGGTTGCACAAGTCTCGGCGATTATCAACTGCAAAAGTTGATGACCAGTACGTTCGATACCGGGAAGGGCCATGAGGGTTTACACATCGAGACTATAACGCGCGCCTAACGACACAAGAAGATA",
                "GTTAGCACAGTCGGATGACTGGTATAATGTTCGTCAACGATTACATAGGACTCGACTGATCGGGACGAACTTAGGCATAAGGGGGAAAGTAGCCTCCCTTCTATACCGGAAGGTTTGGTATACTCTCCCCGTCTGGCAGAATCGTCCCTCTATTTTGTCCACTAGTTTACACGTGGCAGAAGCCGGCCGGGGTATTGCCTCGACTCGCTATGCGAAGTGGAGGCCAAACGATAACTAAATAAACAAGGACCACTACGGATGTGAGTGGGCACCACTTAAGCGTCTGCATCGAACACATGGGAAAACCTTTCGTTAAAGCAAACCTAAAGTTAGACCACGGACCTCTTGTGATGCTAATTCGGACGGAGGGCCAGCCAAGGGAGTCTGGGTTACGGGTCTTGGGTAGACCTTTCAGTCAATGCCTCCATTTTGGTCCAAACAACGAATGTTCAATCACTCGGTAGTAGCCCAGTGTTATGAGTGGACATCAGTTGACCAGGATGCCCACGCCTATGGTAACGGCGAGTGTAGAGTCCGAATAGGTGGTATCCTATCGAGACTGTTTGTACGTAAAACCGGCATTCAGCAACTGTTGACCGGTGGTAATCCTCGAGGGGTCTAATCACCGATCCATTAGAGTACCCTGTACTACCCTATACAACCGAAGGTAAATCGATACCTAACAGGGGTTATGCGCTCGTAAATGCTCCCCGGAGTGTACGCTGCTCGGCTAAGGCGTTCCAAACTTGTTAAAATCTTTACTCGGATTACTTGATGGGACGATTTACCACACCCACAGGCTCATCTATACGTTGATCGGAGCTGCGATTAACAGACGGTAACGCACGGTGGTATCAAAACATTACATTCATTAAGCATGAGGGTGTCACATTGAAAGTACCAGGCATAAAAATTGCCTGTAACCTTGGGATTCGTAAATCCAGCCGAAGGCTGACAAATTCGGGCGGAGGACACTAATTTACAACACCTCGGAGGTAAG"
		};
		
		seqNr=seqNr%sequencesArray.length;
		
		final String seqstr=sequencesArray[seqNr];
		Sequence seq=new Sequence(taxonID, seqstr);
		final TaxonSet txs=tret.m_taxonset.get();
		final Input<TaxonSet> txset=tret.m_taxonset;
		Alignment ali=txs.alignmentInput.get();

		// final int [] distances=new int[nrOfSequencessBeforeUpdate];
		// get the sequence of the last leaf
		final int nrOfElements=seqstr.length();
		
		Arrays.parallelSetAll(distances, e ->
       	{ 
       		String curSeqString=ali.getSequence(e).dataInput.get();
       		int dist=0;
       		// calculate the distance
       		for(int i=0;i<nrOfElements;i++)
       		{
       			if(seqstr.charAt(i)!=curSeqString.charAt(i))
       			{
       				dist++;
       			}
       		}
       		return dist;
       	});
		
		// check if there is any distance equal to 0 (implies that this sequence was already in)
		boolean sequenceAlreadyPresent = Arrays.stream(distances).parallel().anyMatch(i -> i==0);
		
		if(sequenceAlreadyPresent)
		{// signal that the sequence is already in the sequence list
			
			return null;
		}
		
		//List<Sequence> 
		ali.initializeWithAddedSequenceList(Arrays.asList(seq), false);
		//txs.addTaxaName(taxonID);
		final Input<Alignment> aliinput=txs.alignmentInput;
		

		final Integer lengthOfSeq=nrOfElements;
		final double lenSeq=lengthOfSeq.doubleValue();
		
		// needed for drawing of the height
		final double sdOfGaussian=1.0/Math.sqrt(lenSeq);

		
		// after having updated the tree we need to update all obj that hv the tree as input?
		// probably not needed as the initialisation is done anyway in the mcmc init
		
		// the first element is the leaf, the second is the distance of the new leaf from each of the previous ones
		HashMap<Integer, Integer>  selectedLeafMap = new HashMap<Integer, Integer>();
		for(int i=0; i< beastMClist.length; i++)
		{// initialise
			selectedLeafMap.put(i, -1);
		}

		Arrays.parallelSetAll(beastMClist, e ->
	       	{ 
	       		// add sequence here to all
	       		// would be good to have a shared taxa, taxonset, alignment for all particles
	       		
	       		MCMC mc=beastMClist[e].m_mcmc;
	       		State stt=mc.getState();
	       		double theta=stt.stateNode[popsizepositionInStateArray].getArrayValue();

	       		/* start of the part of drawing the leaf */
				Integer selectedLeaf=0;
				//double theta=beastMClist[e];
				{
					double qtToElevate=(lenSeq*theta)/(nrOfSequencessBeforeUpdate+(lenSeq*theta)); // this will hv to change and contain the effective pop size
					double []probabilityWeight= new double[distances.length];//{0.1,0.2,0.25,0.85};//
					Arrays.parallelSetAll(probabilityWeight, ee -> {return Math.exp(distances[ee]*Math.log(qtToElevate));});
					org.apache.commons.math3.util.Pair<Integer, Double> itemToInit=new org.apache.commons.math3.util.Pair<Integer, Double>(0,0.0);
					List<org.apache.commons.math3.util.Pair<Integer, Double>> pmfWeights=new ArrayList<org.apache.commons.math3.util.Pair<Integer, Double>>(Collections.nCopies(probabilityWeight.length, itemToInit));
	        		
					Arrays.parallelSetAll(probabilityWeight, ee ->{
	        			pmfWeights.set(ee, new org.apache.commons.math3.util.Pair<Integer, Double>(ee,probabilityWeight[ee]));
	        			return probabilityWeight[ee];
	        		});
					EnumeratedDistribution enDist=new EnumeratedDistribution<>(pmfWeights);
					// the following is the draw of the leaf (position of the leaf in the array)
					selectedLeaf=(Integer) enDist.sample();
					selectedLeafMap.replace(e, selectedLeaf);
				}
	       		/* end of draw of the leaf */
	       		double my_height=0.0;
				double meanOfGaussian=2*Math.asin(Math.sqrt(distances[selectedLeaf.intValue()]/lenSeq));
	       		/* start of the part of drawing the height */
				{
					double upperBoundTruncatedGaussian=2.094395102393195; // this is 2*Math.asin(sqrt(3)/2.0);
					//org.apache.commons.math3.distribution.NormalDistribution norm=new org.apache.commons.math3.distribution.NormalDistribution(meanOfGaussian, sdOfGaussian);
					//TruncatedNormal tn = new TruncatedNormal(meanOfGaussian, sdOfGaussian,  Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian);

					double inerval;
						double beta=TruncatedNormal.sampleUpgraded(Double.NEGATIVE_INFINITY, upperBoundTruncatedGaussian); //tn.sample();
						if(beta > upperBoundTruncatedGaussian)
						{
				            throw new RuntimeException(
				                    "Unable to draw properly from the truncated Gaussian\n");
						}
						// from the paper on Sequential Monte Carlo transformations,
						// calculate the height (formula 24 of paper)
						// decide if it is better to have a different height for every particle:
						// the process of selectin height is independent of the tree
						double sinVal=Math.sin(beta/2.0);
						inerval=1.0-((4.0/3.0)*sinVal*sinVal);

						if(inerval > 0)
						{
						    my_height=(-3.0/((4.0)*theta))*Math.log(1.0-((4.0/3.0)*sinVal*sinVal));
						}
						else
						{
				            throw new RuntimeException(
				                    "Error in formula from the truncated Gaussian, negative log argument!!!\n");
						}
				}
	       		/* end of draw of the height */

				
				/* 
				 * here we can calculate the components for the weights
				 */
				/*
				 * Leaf first
				 */
				
				// we need to calculate these quantities after the tree has been updated (i.e. after sequence added
				
	       		   double Ms=distances[selectedLeaf.intValue()];
	       		   double N=lenSeq; // length of sequence
	       		   double t=nrOfSequencessBeforeUpdate;
				   double logUnnormalisedIncrementalWeightsLeafSelection=-Ms*(Math.log(N*theta)-Math.log(t+N*theta));
				/*
				 * height afterwards
				 */
				   double logUnnormalisedIncrementalWeightsHeightSelection=calcLogHeightSelection(theta, my_height, meanOfGaussian, sdOfGaussian);
				  
				   
	       		Tree tretre=(Tree) stt.stateNode[treepositionInStateArray];
	       		TaxonSet txs_x=tretre.m_taxonset.get();
	       		
	       		// set same alignment and taxonset to all
	       		// also set startstate for the mcmc

	        	try {
						RandomTree.insertSequenceCoalescent_updated(tretre,selectedLeaf.intValue(), my_height);
			        	List<StateNodeInitialiser> inits=beastMClist[e].m_mcmc.initialisersInput.get();
			        	for(StateNodeInitialiser st:inits)
			        	{
			        		if (st instanceof RandomTree)
			        		{
			        			((RandomTree)st).taxaInput=aliinput;
			        			RandomTree.copyTree((RandomTree)st, tretre);
			        		}
			        	}
			       		txs_x=txs;
			       		tretre.m_traitList=tret.m_traitList;
			       		tretre.m_taxonset=tret.m_taxonset;
			       		tretre.nodeTypeInput=tret.nodeTypeInput;
			       		mc.startStateInput.setValue(stt, null);
			       		CompoundDistribution dst=(CompoundDistribution) mc.GetLikelihoodFromInput();
			       		Map<String, Input<?>> mp=dst.getInputs();
			       		mp.get("Alignment");
			       		List<Distribution> trdt=dst.pDistributions.get();
			       		for(Distribution distrib:trdt)
			       		{
			       			if(distrib instanceof ThreadedTreeLikelihood)
			       			{
				       			ThreadedTreeLikelihood trd=(ThreadedTreeLikelihood)distrib;
				       			trd.dataInput=aliinput;
				       			trd.treeInput.setValue(tretre, null);
			       			}
			       			distrib.initAndValidate();
			       		}
			       		
				} catch (InstantiationException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (IllegalAccessException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (ClassNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} 	       		
	        	return beastMClist[e];
    	}
    	);
		
		return selectedLeafMap;
	}

    private static void addNormalisedWeights(double[] in_logWeightsNormalized, double[] out_logIncrementalWeights)
    {
		Arrays.parallelSetAll(out_logIncrementalWeights, e ->{
			return out_logIncrementalWeights[e]+in_logWeightsNormalized[e];
		});
    }
    
    private static void addSequenceUpdated(Sequential[] beastMClist, int treepositionInStateArray, int popsizepositionInStateArray, int nrOfSequencessBeforeUpdate, final int [] distances, ArrayList<List<Double>> values)
	{
		// add the new sequence first to first element of the list, then use the same shared input for all
		Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
		// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
		// int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
		int numberOfSequencesAfterUpdate=nrOfSequencessBeforeUpdate+1;
		String taxonID="t"+numberOfSequencesAfterUpdate;
		final TaxonSet txs=tret.m_taxonset.get();
		final Input<Alignment> aliinput=txs.alignmentInput;
		
				  
		Arrays.parallelSetAll(beastMClist, e ->
       	{ 
       		// add sequence here to all
       		// would be good to have a shared taxa, taxonset, alignment for all particles
       		
       		MCMC mc=beastMClist[e].m_mcmc;
       		State stt=mc.getState();
	       		Tree tretre=(Tree) stt.stateNode[treepositionInStateArray];
	       		TaxonSet txs_x=tretre.m_taxonset.get();
	       		
	       		// set same alignment and taxonset to all
	       		// also set startstate for the mcmc

	       		try {
						RandomTree.insertSequenceCoalescent_updated(tretre,values.get(e).get(selectedLeafPosition).intValue(), values.get(e).get(heightPosition).doubleValue());
			        	List<StateNodeInitialiser> inits=beastMClist[e].m_mcmc.initialisersInput.get();
			        	for(StateNodeInitialiser st:inits)
			        	{
			        		if (st instanceof RandomTree)
			        		{
			        			((RandomTree)st).taxaInput=aliinput;
			        			RandomTree.copyTree((RandomTree)st, tretre);
			        		}
			        	}
			       		txs_x=txs;
			       		tretre.m_traitList=tret.m_traitList;
			       		tretre.m_taxonset=tret.m_taxonset;
			       		tretre.nodeTypeInput=tret.nodeTypeInput;
			       		mc.startStateInput.setValue(stt, null);
			       		CompoundDistribution dst=(CompoundDistribution) mc.GetLikelihoodFromInput();
			       		Map<String, Input<?>> mp=dst.getInputs();
			       		mp.get("Alignment");
			       		List<Distribution> trdt=dst.pDistributions.get();
			       		for(Distribution distrib:trdt)
			       		{
			       			if(distrib instanceof ThreadedTreeLikelihood)
			       			{
				       			ThreadedTreeLikelihood trd=(ThreadedTreeLikelihood)distrib;
				       			trd.dataInput=aliinput;
				       			trd.treeInput.setValue(tretre, null);
			       			}
			       			distrib.initAndValidate();
			       		}
			       		
				} catch (InstantiationException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (IllegalAccessException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (ClassNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} 	       		
	        	return beastMClist[e];
    	}
    	);
		
	}
    
    public static void main(String[] args) {
 
    	System.out.println("Main 0.....");
    	long sleepseconds;
    	boolean takeExponentsFromFile;
    	String exponentFileName;
    	
        { // nr of SMC particles, we need to know rightaway
        	Arguments arguments=parseArguments(args);
            if (arguments.hasOption("nrOfSMCparticles")) {
            	NR_OF_PARTICLES=arguments.getLongOption("nrOfSMCparticles");
            }
            else
            {
            	NR_OF_PARTICLES=1;
            }
            System.out.println("If it is an SMC run, I will be using "+ NR_OF_PARTICLES+" particles");
        }
        
        { // nr of SMC particles, we need to know rightaway
        	Arguments arguments=parseArguments(args);
            if (arguments.hasOption("sleepseconds")) {
            	sleepseconds=arguments.getLongOption("sleepseconds");
            }
            else
            {
            	sleepseconds=0;
            }
        }

        {
        	Arguments arguments=parseArguments(args);
        	if(arguments.hasOption("targetCESS"))
        	{
        		targetCESSval=arguments.getRealOption("targetCESS");
        	}
        	else
        	{
        		targetCESSval=0.9;
        	}
            System.out.println("Target CESS is: "+ targetCESSval);
        }

        {
        	Arguments arguments=parseArguments(args);
            if (arguments.hasOption("exponentFile")) 
            {
            	takeExponentsFromFile=true;
            	exponentFileName = arguments.getStringOption("exponentFile");
            }
            else
            {
            	takeExponentsFromFile=false;
            	exponentFileName="";
            }
        }
    	
        {
        	Arguments arguments=parseArguments(args);
            if (arguments.hasOption("adaptiveVariance")) 
            {
            	adaptiveVariance=true;
            }
            else
            {
            	adaptiveVariance=false;
            }
        }

        try {
            System.setProperty("beast.debug", "true"); // 
            
            long N=BeastMCMC.NR_OF_PARTICLES;
            
            int N_int= (int) N; // to avoid repeated casting
            
            // constant set to log(1/N) for practical purposes (we don't need to recalculate it over n over)
            final double minuslogN=-java.lang.Math.log(N); // log(1/N) using log properties

            // variable for the annealing, how many steps to arrive from 0 to 1 (for ex. if 100 then steps are 0.01, 0.02...)
            final int maxvalcnt=500; // this is nr of steps minus 1
            
        	// variables for the weights in log space
        	double logIncrementalWeights[] = new double[N_int]; // vector of weights for the particles
            double logWeightsNormalized[] = new double[N_int]; // vector of weights for the particles
            
            //BeastDialog dlg=CreateAndShowDialog();
            //double currentExponent, previousExponent; //exponent to be used for simulated annhealing
            long exponentCnt;
            Path currentRelativePath = Paths.get("");
            String s = currentRelativePath.toAbsolutePath().toString();

            // init the random generator: MUST do otherwise each run might bring same results
            //initRandomizer();

        	// this is the list of the particles of the SMC
            Sequential[] beastMClist = new Sequential[N_int];
            if(sleepseconds>0)
            {
            	Thread.sleep(sleepseconds*1000);
            }

           	// here we init the list of particles and sample from the prior
        	initParticlesAndSampleFromPrior(beastMClist, null, args);
            
            // get the position of the tree in the state array
            int treepositionInStateArray=getTreePosition((MCMC)beastMClist[0].m_runnable); // position of the tree in the state vector

            if(treepositionInStateArray<0)
            {
            	System.err.println("Unable to find Tree class in statenode");
                System.exit(0);
            }

            // get the position of the tree in the state array
            int populationsizePositionInStateArray=getPopsizePosition((MCMC)beastMClist[0].m_runnable); // position of the tree in the state vector

            if(populationsizePositionInStateArray<0)
            {
            	System.err.println("Unable to find Tree class in statenode");
                System.exit(0);
            }
            
           	Arrays.parallelSetAll(logWeightsNormalized, e->{
           		Sequential bmcc=beastMClist[e];
           		bmcc.m_mcmc.setPopSizePositionInStateArray(populationsizePositionInStateArray);
           		return logWeightsNormalized[e];
           	});
            
            
            // save the ESS
            ByteArrayOutputStream ess = new ByteArrayOutputStream();
            PrintStream outEss = new PrintStream(ess);

            // save the CESS
            ByteArrayOutputStream cess = new ByteArrayOutputStream();
            PrintStream outCEss = new PrintStream(cess);

            // save the normalised weights
            ByteArrayOutputStream weightsStream = new ByteArrayOutputStream();
            PrintStream outWeights = new PrintStream(weightsStream);
            
           	// in the following initialize the weights to 1/N
            initNotmalisedWeights(logWeightsNormalized, minuslogN);
			
        	//currentExponent=0.0;            	
        	double ESSval, CESSval;
        	
        	// string used within the files as row counter
        	String rowCounterString;
        	final String divider=",";
        	long startTime=System.nanoTime();
        	long[] nrOfMCMCrejections=new long[N_int];
        	//List<Integer> nrOfMCMCrej=new ArrayList<>();
        	List<Double> avgRejectionList=new ArrayList<>();
        	double[] avgRejection=new double[maxvalcnt];
        	
        	// init all the avgRejection elements to "not done"
        	final long MCMC_NotDone=-1;
            Arrays.parallelSetAll(avgRejection, e->MCMC_NotDone);
        	double auxDoubleVar;
            
        	// size of the step for the exponent
        	double stepSize=1/((double)maxvalcnt);
        	
        	/*
        	 static int getCESSexponent(Sequential[] beastMClist, double[] logIncrementalWeights, double[] logWeightsNormalized, double stepSize, double previousExponent, int exponentCnt, int maxvalcnt, double outNextExponent)
 
        	*/
        	
        	// here add the new sequence
        	
        	boolean useCESS=false;
        	int nextCESS;
            
        	double nextExponentDouble=0.0,currentExponentDouble=0.0;
        	
        	// list of exponents, to be used if there is a input argument to take predefined exponents
        	List<Double>expList=null;
        	
        	if(takeExponentsFromFile == true)
        	{
        		  expList= new ArrayList<>();
        		  //String exponentFileName="/Users/lr411/Leo/Github/xmlFilesBeast/CESS_20200103_030123_P500_E218_S100000.txt";
        		  try
        		  {
        		    BufferedReader reader = new BufferedReader(new FileReader(exponentFileName));
        		    String line;
        		    while ((line = reader.readLine()) != null)
        		    {
        		      String[] rcrds=line.split(",");
        		      if(rcrds[0] != null)
        		    	  expList.add(Double.parseDouble(rcrds[0]));
        		    }
        		    reader.close();
        		  }
        		  catch (Exception e)
        		  {
        		    System.err.format("Exception occurred trying to read '%s'.", exponentFileName);
        		    e.printStackTrace();
        		  }
        	}
        	
        	
        	// this variable is used only if u take the exponent from a file
            int exponentIndex=0;
            
            // this variable is used only if you want to add sequences on the fly
            boolean changeSeq=false;
            
            while(nextExponentDouble<1)//for (exponentCnt=0; exponentCnt<maxvalcnt; exponentCnt++)
            {// starts from the prior and goes to target (reached when the exponent is equal to 1)
            	// smcStates[(int)i][(int)exponentCnt]=mc.getState();
            	//if(exponentCnt>=(maxvalcnt/100))
            	//	break;
            	
            	double outNextExponent=0.0;
//           	    int nextCESS=getCESSexponent(beastMClist, logIncrementalWeights, logWeightsNormalized, stepSize, previousExponent, (int) exponentCnt, maxvalcnt+1, outNextExponent);
                if(sleepseconds>0)
                {
                	Thread.sleep(sleepseconds*1000);
                }
                
				if(useCESS)
				{
				   if(takeExponentsFromFile)
				   {
					   if(exponentIndex>expList.size())
					   {
						   nextExponentDouble=1.0;
					   }
					   else
					   {
						   nextExponentDouble=expList.get(exponentIndex).doubleValue();
						   exponentIndex++;
					   }
				   }
				   else
				   {
					   nextExponentDouble=getCESSexponent_double(beastMClist, logIncrementalWeights, logWeightsNormalized, currentExponentDouble, targetCESSval);
					   if(nextExponentDouble>1.0)
						   nextExponentDouble=1.0;
				   }
				}
				else
				{
					nextExponentDouble=nextExponentDouble+stepSize;
				}
				
				if(currentExponentDouble==0)
				{
					int dbgvar=0;
				}
				
				System.out.println("Exponent: "+nextExponentDouble);
				
                if(sleepseconds>0)
                {
                	Thread.sleep(sleepseconds*1000);
                	sleepseconds=0;
                }
				// reweight done below, calculation of the incremental part
            	calculateIncrementalWeights(beastMClist, logIncrementalWeights, currentExponentDouble, nextExponentDouble);
				// CESS to be calculated before renormalising
            	
            	currentExponentDouble= nextExponentDouble;

            	CESSval=CESS(logIncrementalWeights, logWeightsNormalized);

				System.out.println("CESS: "+CESSval);

				rowCounterString=currentExponentDouble + divider;
            	
				outCEss.println(rowCounterString + CESSval);

				// double cessNorm=CESS_Normalised(logIncrementalWeights, logWeightsNormalized);
                // normalising below, logWeightsNormalized is the output, logIncrementalWeights the input
               	normaliseWeights(logIncrementalWeights, logWeightsNormalized);
               	
               	if(adaptiveVariance)
               	{
                   	// here calculate current mean and variance of the set of particles for populationSize, we use the weights
                   	final double currentVar=calculateParameterVariance(beastMClist, logWeightsNormalized, populationsizePositionInStateArray);
                   	//setParticleVariance
                   	Arrays.parallelSetAll(logWeightsNormalized, e->{
            			//BeastMCMC bmc=beastMClist[e];
                   		Sequential bmcc=beastMClist[e];
                   		bmcc.m_mcmc.setParticleVariance(currentVar);
                   		return logWeightsNormalized[e];
                   	});

               	}

				List <Integer> times = new ArrayList<>();
				
				if(changeSeq)
				{
					if(nextExponentDouble>0.05)
					{
						changeSeq=false;
			        	long starttimeNano;
						long elapsedTimeNano;
			        	long elapsedTimeMs;
			        	boolean calledBeforeTreeUpdated;

						for(int loc=0;loc<10;loc++)
						{
							Tree tret=(Tree)beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray];
							// the numberOfSequencesAfterUpdate is return value: number of updated number of leaves
							State st=beastMClist[0].m_mcmc.getState();
							int nrOfSequencessBeforeUpdate=tret.getLeafNodeCount();
							// the distances array will be calculated inside the routine to add the sequence
							int [] distances=new int[nrOfSequencessBeforeUpdate];
							
							// check if the sequence is already present
						    if(checkAndAddSequenceAndCalculateDistances(beastMClist, treepositionInStateArray, m_sequencesArray[loc%m_sequencesArray.length], nrOfSequencessBeforeUpdate, distances))
						    {// all good, the sequence was not already in the data, proceed

						        int sequenceLength=m_sequencesArray[0].length();
						        ArrayList<HashSet<Integer>> leaves = new ArrayList<>();
								for(int i=0; i< beastMClist.length; i++)
								{// initialise and add a list of doubles for
									leaves.add(new HashSet<Integer>());
								}
						    	ArrayList<List<Double>> retvals=drawLeafAndHeightAndCalculateLogUpdated( beastMClist, treepositionInStateArray, populationsizePositionInStateArray,  sequenceLength, nrOfSequencessBeforeUpdate,  distances, leaves);
						    	calledBeforeTreeUpdated=true;
							    calculateIncrementalWeightsForTransformation(beastMClist, logIncrementalWeights, currentExponentDouble, calledBeforeTreeUpdated);
								starttimeNano=System.nanoTime();
								// add the new sequence first to first element of the list, then use the same shared input for all
								
								//HashMap<Integer,Integer> selectedLeaves=addSequence(beastMClist, treepositionInStateArray, populationsizePositionInStateArray, loc, nrOfSequencessBeforeUpdate, distances);			        		
								addSequenceUpdated(beastMClist, treepositionInStateArray, populationsizePositionInStateArray, nrOfSequencessBeforeUpdate, distances, retvals);			        		
								
								elapsedTimeNano=System.nanoTime()-starttimeNano;
								//elapsedTimeMs=getExecutionLength(elapsedTimeNano);
								times.add((int)(elapsedTimeNano/1000000));
								// here update transformation weight
								calledBeforeTreeUpdated=false;
							    calculateIncrementalWeightsForTransformation(beastMClist, logIncrementalWeights, currentExponentDouble, calledBeforeTreeUpdated);
							    calculateIncrementalWeightsForTransformationAddLeafComponent(beastMClist, logIncrementalWeights, retvals);
							    calculateIncrementalWeightsForTransformationAddHeightComponent(beastMClist, logIncrementalWeights, retvals);
							    addNormalisedWeights(logIncrementalWeights, logWeightsNormalized);
							    // here normalise weights
				               	normaliseWeights(logIncrementalWeights, logWeightsNormalized);
						    }
						}
					}
				}

				String strng=Arrays.toString(logIncrementalWeights);
               	outWeights.println(strng);
               	
               	ESSval=ESS(logWeightsNormalized);
				
               	System.out.println("ESS: "+ESSval);
               	
				outEss.println(rowCounterString + ESSval);
				

               	if(nextExponentDouble<=1.0) 
               	{
    				if((ESSval<(N_int/2.0)) || (nextExponentDouble==1.0))
                   	{// at the moment we resample always at the end to make it like the mcmc
        				List<Integer> stratifiedList=stratified_resample((double [])logWeightsNormalized);

                   		// update the particles list according to the list output from the resample process
                   		updateParticlesList(stratifiedList, beastMClist);
                       	// in the following initialize the weights to 1/N
                        initNotmalisedWeights(logWeightsNormalized, minuslogN);  
                   	}

                   	// only do MCMC move if we are not in the last step
                   	//if(exponentCnt<maxvalcnt-1) 
                   	{
                        // do the mcmc moves on the particles and set the exponent for annealing
                        doMCMC_andSetExponentForAnnealing(beastMClist, currentExponentDouble, nrOfMCMCrejections);
                        auxDoubleVar=Arrays.stream(nrOfMCMCrejections).average().getAsDouble();
                        avgRejectionList.add(auxDoubleVar);
                        //avgRejection[(int) exponentCnt]=auxDoubleVar;
                   	}
/*    				System.out.println("After");
    				for(int i=0; i<beastMClist.length; i++)
    				{
    					Tree tr=((Tree) beastMClist[0].m_mcmc.getState().stateNode[treepositionInStateArray]);
    					System.out.println("particle "+i+", nodes: "+tr.getNodeCount()+", stored: "+tr.getStoredNodes().length);
    				}
*/
               	}
            }// outer parentheses 
        	
        	// how much time did it take
        	long elapsedTimeNano=System.nanoTime()-startTime;
        	long elapsedTime=getExecutionLength(elapsedTimeNano);
        	
        	// put information 
            String informativeAppendString=formatAppendString(N_int, elapsedTime, maxvalcnt);

            // save the logs with parameters
            saveLogs(informativeAppendString, ess, cess, weightsStream);
            
            // saveStateSpaceLogs(informativeAppendString,popSizeStream, gammaShapeStream);
            
            // save average rejections
            saveAvgRej(informativeAppendString, avgRejectionList, "AvgRejection");
            
            // save the tree particles
            //saveTreeParticles(beastMClist, informativeAppendString, treepositionInStateArray, N_int);
            
            // save state space array
            saveStateSpaceParticles(beastMClist, informativeAppendString, N_int);
            
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(BeastMCMC.getUsage());
        } catch (Throwable e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
        if (System.getProperty("beast.useWindow") == null) {
            // this indicates no window is open
            System.exit(0);
        }
    } // main


} // BeastMCMC
