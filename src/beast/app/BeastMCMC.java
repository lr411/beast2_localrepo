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
import beast.app.beastapp.BeastDialog;
import beast.app.beastapp.BeastMain;
import beast.app.beauti.Beauti;
import beast.app.draw.ExtensionFileFilter;
import beast.app.util.Version;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Runnable;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.util.*;
import jam.util.IconUtils;
import org.json.JSONException;
import org.xml.sax.SAXException;

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
    public static final long NR_OF_PARTICLES = 500;
    // path to save the logs
    final static String logsPath="/Users/lr411/Leo/Github/Genomics/logs_BEAST2/";
    // nr of MCMC moves
    public static final int NR_OF_MCMC_MOVES = 5;
    public BEASTInterface m_treeobj=null;
    public double m_initPopSize;
    public double m_gammaShapeLog;
    public double m_gammaShape;
    /**
     * number of threads used to run the likelihood beast.core *
     */
    static public int m_nThreads = 1;
    
    //static private boolean m_dialogInitialized=false;
    /**
     * thread pool *
     */
    public static ExecutorService g_exec = Executors.newFixedThreadPool(m_nThreads);
    /**
     * random number seed used to initialise Randomizer *
     */
    long m_nSeed = 127;

    BeastDialog m_dialog=null;
    
    /**
     * MCMC object to execute *
     */
    Runnable m_runnable;

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
    public static void initRandomizer(long seed)
    {    	
        Randomizer.setSeed(seed);
    }
    
    // init the randomizer using the current time in milliseconds
    public static void initRandomizer()
    {
	    Calendar calendar = Calendar.getInstance();
	    long randSeed = calendar.getTimeInMillis();
    	
        Randomizer.setSeed(randSeed);
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
        final Version version = new BEASTVersion2();
        final String titleString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                "Version " + version.getVersionString() + ", " + version.getDateString() + "</p></center></html>";
        final javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");
        final String nameString = "BEAST " + version.getVersionString();
        
        final BeastDialog dialog = new BeastDialog(new JFrame(), titleString, icon);
        
       // if(m_dialogInitialized == false)
        {
	        if (!dialog.showDialog(nameString, 1)) {
	            return null;
	        }
        }
        
       // m_dialogInitialized=true;

        return dialog;
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
                            throw new IllegalArgumentException("Wrong argument");
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error parsing command line arguments: " + Arrays.toString(args) + "\nArguments ignored\n\n" + getUsage());
        }

        if (beastFile == null) {
            // Not resuming so get starting options...

            
           // boolean oldDialogInitialised = this.m_dialogInitialized;
            
            
          //  if(oldDialogInitialised == false)
        //    {// only do once
                // Leo: create and show only once
                //m_dialog = CreateAndShowDialog();
        	    BeastDialog dialog=m_dialog; // just to avoid changing names

                List<String> MCMCargs = new ArrayList<>();

                switch (m_dialog.getLogginMode()) {
	                case 0:/* do not ovewrite */
	                    break;
	                case 1:
	                    MCMCargs.add("-overwrite");
	                    break;
	                case 2:
	                    MCMCargs.add("-resume");
	                    break;
	            }
	            MCMCargs.add("-seed");
	            MCMCargs.add(m_dialog.getSeed() + "");
	
	            if (m_dialog.getThreadPoolSize() > 0) {
	                MCMCargs.add("-threads");
	                MCMCargs.add(m_dialog.getThreadPoolSize() + "");
	            }
          //  }
	            
            boolean useBeagle = dialog.useBeagle();
            boolean beagleShowInfo = false;
            long beagleFlags = 0;
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
            if (beagleFlags != 0) {
                System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));
            }
            if (!useBeagle) {
                System.setProperty("java.only", "true");
            }

            File inputFile = dialog.getInputFile();
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

    public static double ESS(double[] normalized_log_weights)
    {
    	double[] weightsSquared= new double[normalized_log_weights.length];//Arrays.copyOf(normalized_log_weights, normalized_log_weights.length);
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
	           		Class cls=stn.getClass();
	           		if(cls.getName().endsWith(".Tree"))
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
       	    MCMC mc=bmcc.m_mcmc;

       	    // the mcmc run is done with the previous exponent
        	try {
				bmcc.run();
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

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
    
    // the output is the updated incremental weights
    public static void calculateIncrementalWeights(Sequential[] beastMClist, double[] output_logUnnormalisedIncrementalWeights, final double previousExponent, final double nextExponent)
    {
       	Arrays.parallelSetAll(output_logUnnormalisedIncrementalWeights, e->{
			//BeastMCMC bmc=beastMClist[e];
			MCMC mc=beastMClist[e].m_mcmc;
			// performance wise, we don't need the log1 and log 2 we could substitute the full expressions to log1 and 2
			double log1=mc.calculateLogPSimulatedAnnhealing(nextExponent);
			double log2=mc.calculateLogPSimulatedAnnhealing(previousExponent);
//        	// ?? is it ok to set the exponent here????
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
            bmc.SetDlg(dlg);
       	    // the mcmc run is done with the previous exponent
        	try {
                bmc.parseArgs(args);
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
            MCMC mc= (MCMC)bmc.m_runnable;
            mc.SetDistributionsFromInput();
            // initialise the state of the posterior
            mc.initStateAndPosterior();
            
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
    
    public static void main(String[] args) {
        
    	try {
            System.setProperty("beast.debug", "true"); // 
            
            long N=BeastMCMC.NR_OF_PARTICLES;
            
            int N_int= (int) N; // to avoid repeated casting
            
            // constant set to log(1/N) for practical purposes (we don't need to recalculate it over n over)
            final double minuslogN=-java.lang.Math.log(N); // log(1/N) using log properties

            // variable for the annealing, how many steps to arrive from 0 to 1 (for ex. if 100 then steps are 0.01, 0.02...)
            final int maxvalcnt=100; // this is nr of steps minus 1
            
        	// variables for the weights in log space
        	double logIncrementalWeights[] = new double[N_int]; // vector of weights for the particles
            double logWeightsNormalized[] = new double[N_int]; // vector of weights for the particles
            
            BeastDialog dlg=CreateAndShowDialog();
            double currentExponent, previousExponent; //exponent to be used for simulated annhealing
            long exponentCnt;
            Path currentRelativePath = Paths.get("");
            String s = currentRelativePath.toAbsolutePath().toString();

            // init the random generator: MUST do otherwise each run might bring same results
            initRandomizer();

        	// this is the list of the particles of the SMC
            Sequential[] beastMClist = new Sequential[N_int];

           	// here we init the list of particles and sample from the prior
        	initParticlesAndSampleFromPrior(beastMClist, dlg, args);
            
            // get the position of the tree in the state array
            int treepositionInStateArray=getTreePosition((MCMC)beastMClist[0].m_runnable); // position of the tree in the state vector

            if(treepositionInStateArray<0)
            {
            	System.err.println("Unable to find Tree class in statenode");
                System.exit(0);
            }
            
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
			
        	currentExponent=0.0;            	
        	double ESSval, CESSval;
        	
        	// string used within the files as row counter
        	String rowCounterString;
        	final String divider=",";
        	long startTime=System.nanoTime();
        	long[] nrOfMCMCrejections=new long[N_int];
        	double[] avgRejection=new double[maxvalcnt];
        	
        	// init all the avgRejection elements to "not done"
        	final long MCMC_NotDone=-1;
            Arrays.parallelSetAll(avgRejection, e->MCMC_NotDone);
        	double auxDoubleVar;
            
            for (exponentCnt=0; exponentCnt<maxvalcnt; exponentCnt++)
            {// starts from the prior and goes to target (reached when the exponent is equal to 1)
            	// smcStates[(int)i][(int)exponentCnt]=mc.getState();
            	previousExponent=currentExponent;                	
            	
            	currentExponent=((double)exponentCnt+1)/((double)maxvalcnt);
            	// updates the weights with the ratio of future-current functions calculated at the current state
				
            	rowCounterString=exponentCnt + divider + currentExponent + divider;
            	
               	// reweight done below, calculation of the incremental part
            	calculateIncrementalWeights(beastMClist, logIncrementalWeights, previousExponent, currentExponent);
               	
				// CESS to be calculated before renormalising
            	CESSval=CESS(logIncrementalWeights, logWeightsNormalized);
				outCEss.println(rowCounterString + CESSval);

                // normalising below, logWeightsNormalized is the output, logIncrementalWeights the input
               	normaliseWeights(logIncrementalWeights, logWeightsNormalized);
                String strng=Arrays.toString(logIncrementalWeights);
               	outWeights.println(strng);
               	
               	ESSval=ESS(logWeightsNormalized);
				
				outEss.println(rowCounterString + ESSval);
				
				// only resample if ESS<half particles and if we are not in the last step
               	if(ESSval<(N_int/2.0) && exponentCnt<maxvalcnt-1) 
               	{
    				List<Integer> stratifiedList=stratified_resample((double [])logWeightsNormalized);

               		// update the particles list according to the list output from the resample process
               		updateParticlesList(stratifiedList, beastMClist);
                   	// in the following initialize the weights to 1/N
                    initNotmalisedWeights(logWeightsNormalized, minuslogN);
    				
                    // do the mcmc moves on the particles and set the exponent for annealing
                    doMCMC_andSetExponentForAnnealing(beastMClist, currentExponent, nrOfMCMCrejections);
                    auxDoubleVar=Arrays.stream(nrOfMCMCrejections).average().getAsDouble();
                    avgRejection[(int) exponentCnt]=auxDoubleVar;
               	}
            }// outer parentheses 
        	
        	// how much time did it take
        	long elapsedTimeNano=System.nanoTime()-startTime;
        	long elapsedTime=getExecutionLength(elapsedTimeNano);
        	
        	// put information 
            String informativeAppendString=formatAppendString(N_int, elapsedTime, maxvalcnt);

            // save the logs with parameters
            saveLogs(informativeAppendString, ess, cess, weightsStream);
            
            // save average rejections
            saveParameter(informativeAppendString, avgRejection, "AvgRejection");
            
            // save the tree particles
            saveTreeParticles(beastMClist, informativeAppendString, treepositionInStateArray, N_int);
            
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
