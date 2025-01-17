"""ldmxcfg

Basic python configuration for ldmx-sw application
"""

from pathlib import Path


class EventProcessor:
    """An EventProcessor object

    This object contains the parameters that are necessary for a framework::EventProcessor to be configured.

    You should NOT use this class directly. Use one of the derived classes Producer or Analyzer for clarity.

    Parameters
    ----------
    instanceName : str
        Name of this copy of the producer object
    className : str
        Name (including namespace) of the C++ class that this processor should be
    moduleName : str
        Name of module the C++ class is in (e.g. Ecal or SimCore)
        or full path to the library that should be loaded

    Attributes
    ----------
    histograms : list of histogram1D objects
        List of histogram configure objects for the HistogramPool to make for this processor

    See Also
    --------
    LDMX.Framework.ldmxcfg.Producer : Producer configuration object
    LDMX.Framework.ldmxcfg.Analyzer : Analyzer configuration object
    LDMX.Framework.histogram.histogram : histogram configuration object
    """

    def __init__(self, instanceName, className, moduleName):
        self.instanceName=instanceName
        self.className=className
        self.histograms=[]

        if moduleName.endswith('.so'):
            # assume user passed full path to library
            Process.addLibrary(moduleName)
        else:
            # assume user passed name of module processor is compiled into
            Process.addModule(moduleName)


    @classmethod
    def from_file(cls, source_file, class_name = None, needs = [], instance_name = None, compile_notice = True, **config_kwargs):
        """Construct an event processor "in place" from the passed source file

        Since Framework dynamically loads libraries containing processors after
        the python script has been fully run, we can compile a single-file processor
        into its own library that can then be loaded and run. This function puts
        the library next to the source file and only re-compiles if the source file's
        last modified time is newer than the library (or the library does not exist).

        Note
        ----
        Developing processors in this way is incredible inefficient, especially since it
        does not allow for code to be well organized and split across many files nor does it
        allow for two processors to share common code.
        If you find yourself defining more than one `class` within your new C++ processor,
        it is highly recommended to transition your workflow to including your processor as a
        part of ldmx-sw so that it can fully benefit from a build system.

        Parameters
        ----------
        source_file: str | Path
            path to source file to build into a processor (can be relative to where config is being run)
        class_name: str, default is name of source file
            name of C++ class that is the processor
            defaults to the name of the source file without an extension
        needs: list[str]
            Names of libraries that should be linked to the compiled processor in addition to 'Framework'
            which is linked be default.
            For example, one can gain access to the detector ID infrastructure with 'DetDescr'.
        instance_name: str, default is class_name
            name to give to instance of this C++ processor
        compile_notice: bool, default is True
            print a notice when compilation is triggered
        config_kwargs: dict[str, Any]
            configuration parameters to give to the processor

        Examples
        --------
        A basic walkthrough is available online. https://ldmx-software.github.io/analysis/ldmx-sw.html
        
        If `MyAnalyzer.cxx` contains the class `MyAnalyzer`, then we can put

            p.sequence = [ ldmxcfg.Analyzer.from_file('MyAnalyzer.cxx') ]

        In our config script to run the analyzer on its own in the sequence.
        This default configuration only links the Framework library and so the analyzer
        would only be able to access the Framework and event objects.
        If you needed another library (for example, the 'DetDescr' library has the ID classes),
        one can also

            p.sequence = [ ldmxcfg.Analyzer.from_file('MyAnalyzer.cxx', needs = ['DetDescr']) ]

        To inform the compiler that it should link your analyzer with the 'DetDescr' library.
        **No removal of the library is done** so if you change the `needs` or some other parameter
        to `from_file` you should also remove the library file (`*.so`) before attempting to re-run.

        Returns
        -------
        EventProcessor
            built from the C++ source file and configured with the passed arguments
        """

        if not isinstance(source_file, Path):
            source_file = Path(source_file)
        if not source_file.is_file():
            raise ValueError(f'{source_file} is not accessible.')

        src = source_file.resolve()

        if class_name is None:
            # assume class name is name of file (no extension) if not provided
            class_name = src.stem

        if instance_name is None:
            # use class name for instance name if not provided
            instance_name = class_name

        lib = src.parent / f'lib{src.stem}.so'
        if not lib.is_file() or src.stat().st_mtime > lib.stat().st_mtime:
            if compile_notice:
                print(
                    f'Processor source file {src} is newer than its compiled library {lib}'
                    ' (or library does not exist), recompiling...'
                )
            import subprocess
            libs_to_link = set(['Framework']+needs)
            subprocess.run([
                'g++', '-fPIC', '-shared', # construct a shared library for dynamic loading
                '-o', str(lib), str(src), # define output file and input source file
            ]+[
                f'-l{lib}' for lib in libs_to_link
            ]+[
                '-I/usr/local/include/root', # include ROOT's non-system headers
                '-I@CMAKE_INSTALL_PREFIX@/include', # include ldmx-sw headers (if non-system)
                '-L@CMAKE_INSTALL_PREFIX@/lib', # include ldmx-sw libs (if non-system)
            ], check=True)
            if compile_notice:
                print(f'done compiling {src}')

        instance = cls(instance_name, class_name, str(lib))
        for cfg_name, cfg_val in config_kwargs.items():
            setattr(instance, cfg_name, cfg_val)
        return instance


    def build1DHistogram(self, name, xlabel, bins, xmin = None, xmax = None):
        """Make a 1D histogram 

        If xmin and xmax are not provided, bins is assumed to be
        the bin edges on the x-axis. If they are both provided,
        bins is assumed to be the number of bins on the x-axis.

        Parameters
        ----------
        name : str
            variable name of histogram
        xlabel : str
            title of x-axis of histogram
        bins : int OR list of floats
            Number of bins on x-axis OR bin edges on x-axis
        xmin : float
            Minimum edge of bins on x-axis
        xmax : float
            Maximum edge of bins on x-axis

        See Also
        --------
        LDMX.Framework.histogram.histogram : histogram configuration object
        """

        import LDMX.Framework.histogram as h
        theBinEdges = bins
        if xmin is not None and xmax is not None :
            theBinEdges = h.uniform_binning(bins,xmin,xmax)

        self.histograms.append(h.histogram(name, xlabel,theBinEdges))

    def build2DHistogram(self, name, 
            xlabel = 'X Axis', xbins = 1, xmin = None, xmax = None, 
            ylabel = 'Y Axis', ybins = 1, ymin = None, ymax = None) :
        """Create a 2D histogram

        If {x,y}min or {x,y}max are not provided, {x,y}bins is assumed
        to be the bin edges on the {x,y}-axis. If they are both provided,
        {x,y}-bins is assumed to be the number of bins on the {x,y}-axis.

        Parameters
        ----------
        name : str
            variable name of histogram
        xlabel : str
            title of x-axis of histogram
        xbins : int OR list of floats
            Number of bins on x-axis OR list of bin edges on x-axis
        xmin : float
            Minimum edge of bins on x-axis
        xmax : float
            Maximum edge of bins on x-axis
        ylabel : str
            title of y-axis of histogram
        ybins : int OR list of floats
            Number of bins on y-axis OR list of bin edges on y-axis
        ymin : float
            Minimum edge of bins on y-axis
        ymay : float
            Mayimum edge of bins on y-axis

        See Also
        --------
        LDMX.Framework.histogram.histogram : histogram configuration object

        Examples
        --------

        When doing all uniform binning, you can specify the arguments by position.
            myProcessor.build2DHistogram( 'dummy' ,
                'My X Axis' , 20 , 0. , 1. ,
                'My Y Axis' , 60 , 0. , 10. )

        When using variable binning, you have to use the parameter names.
            myProcessor.build2DHistogram( 'dummy2' ,
                xlabel='My X Axis', xbins=[0.,1.,2.],
                ylabel='My Y Axis', ybins=60, ymin=0., ymax=10. )
        """

        import LDMX.Framework.histogram as h
        theBinEdgesX = xbins
        if xmin is not None and xmax is not None :
            theBinEdgesX = h.uniform_binning(xbins,xmin,xmax)

        theBinEdgesY = ybins
        if ymin is not None and ymax is not None :
            theBinEdgesY = h.uniform_binning(ybins,ymin,ymax)

        self.histograms.append(
                h.histogram(name, xlabel,theBinEdgesX, ylabel,theBinEdgesY)
                )

class Producer(EventProcessor):
    """A producer object.

    This object contains the parameters that are necessary for a framework::Producer to be configured.

    See Also
    --------
    LDMX.Framwork.ldmxcfg.EventProcessor : base class
    """

    def __init__(self, instanceName, className, moduleName):
        super().__init__(instanceName,className, moduleName)

    def __str__(self) :
        """Stringify this Producer, creates a message with all the internal parameters.

        Returns
        -------
        str
            A message with all the parameters and member variables in a human readable format
        """

        msg = "\n  Producer(%s of class %s)"%(self.instanceName,self.className)
        if len(self.__dict__)>0:
            msg += "\n   Parameters:"
            for k, v in self.__dict__.items():
                msg += "\n    " + str(k) + " : " + str(v)

        return msg

class Analyzer(EventProcessor):
    """A analyzer object.

    This object contains the parameters that are necessary for a framework::Analyzer to be configured.

    See Also
    --------
    LDMX.Framework.ldmxcfg.EventProcessor : base class
    """

    def __init__(self, instanceName, className, moduleName):
        super().__init__(instanceName,className, moduleName)

    def __str__(self) :
        """Stringify this Analyzer, creates a message with all the internal parameters.

        Returns
        -------
        str
            A message with all the parameters and member variables in a human readable format
        """

        msg = "\n  Analyzer(%s of class %s)"%(self.instanceName,self.className)
        if len(self.__dict__)>0:
            msg += "\n   Parameters:"
            for k, v in self.__dict__.items():
                msg += "\n    " + str(k) + " : " + str(v)

        return msg

class ConditionsObjectProvider:
    """A ConditionsObjectProvider

    This object contains the parameters that are necessary for a framework::ConditionsObjectProvider to be configured.

    In this constructor we also do two helpful processes.
    1. We append the module that this provider is in to the list of libraries to load
    2. We declare this provider so that the Process "knows" it exists and will load it into the run

    Parameters
    ----------
    objectName : str
        Name of the object this provider provides
    className : str
        Name (including namespace) of the C++ class of the provider
    moduleName : str
        Name of module that this COP is compiled into (e.g. Ecal or EventProc)

    Attributes
    ----------
    tagName : str
        Tag which identifies the generation of information
    """

    def __init__(self, objectName, className, moduleName):
        self.objectName=objectName
        self.className=className
        self.tagName=''

        # make sure process loads this library if it hasn't yet
        Process.addModule(moduleName)
        
        #register this conditions object provider with the process
        Process.declareConditionsObjectProvider(self)

    def setTag(self,newtag) :
        """Set the tag generation of the Conditions

        Parameters
        ----------
        newtag : str
            Tag for generation of conditions
        """

        self.tagName=newtag

    def __eq__(self,other) :
        """Check if two COPs are the same

        We decide that two COPs are 'equal' if they have the same instance and class names
        
        Parameters
        ----------
        other : ConditionsObjectProvider
            other COP to compare agains
        """

        if not isinstance(other,ConditionsObjectProvider) :
            return NotImplemented

        return (self.objectName == other.objectName and self.className == other.className)

    def __str__(self) :
        """Stringify this ConditionsObjectProvider, creates a message with all the internal parameters.

        Returns
        -------
        str
            A message with all the parameters and member variables in a human readable format
        """

        msg = "\n  ConditionsObjectProvider(%s of class %s, tag='%s')"%(self.objectName,self.className,self.tagName)
        if len(self.__dict__)>0:
            msg += "\n   Parameters:"
            for k, v in self.__dict__.items():
                msg += "\n    " + str(k) + " : " + str(v)

        return msg

class RandomNumberSeedService(ConditionsObjectProvider):
    """The random number seed service

    This object registers the random number seed service with the process and
    gives some helper functions for configuration.

    Attributes
    ----------
    seedMode : str
        Name of mode of getting random seeds
    """

    def __init__(self) :
        super().__init__('RandomNumberSeedService','framework::RandomNumberSeedService','Framework')
        self.seedMode = ''
        self.seed=-1 #only used in external mode

        # use run seed mode by default
        self.run()

    def run(self) :
        """Base random number seeds off of the run number"""
        self.seedMode = 'run'

    def external(self,seed) :
        """Input the master random number seed

        Parameters
        ----------
        seed : int
            Integer to use as master random number seed
        """
        self.seedMode = 'external'
        self.seed = seed

    def time(self) :
        """Set master random seed based off of time"""
        self.seedMode = 'time'
    
class Process:
    """Process configuration object

    The python object that stores the necessary parameters for configuring
    a Process for ldmx-app to execute.

    Upon construction, the class-wide reference lastProcess is set
    and the rest of the attributes are set to sensible defaults.

    Parameters
    ----------
    passName : str
        Short reference name for this run of the process

    Attributes
    ----------
    lastProcess : Process
        Class-wide reference to the last Process object to be constructed
    maxEvents : int
        Maximum number events to process.
        If totalEvents is set, this will be ignored.
    maxTriesPerEvent : int
        Maximum number of attempts to make in a row before giving up on an event
        Only used in Production Mode (no input files)
        If totalEvents is set, this will be ignored.
    totalEvents : int
        Number of events we'd like to produce independetly of the number of tries it would take.
        Both maxEvents and maxTriesPerEvent will be ignored. Be warned about infinite loops!
    run : int
        Run number for this process
    inputFiles : list of strings
        Input files to read in event data from and process
    outputFiles : list of strings
        Output files to write out event data to after processing
    sequence : list of Producers and Analyzers
        List of event processors to pass the event bus objects to
    keep : list of strings
        List of rules to keep or drop objects from the event bus
    libraries : list of strings
        List of libraries to load before attempting to build any processors
    skimDefaultIsKeep : bool
        Flag to say whether to process should by default keep the event or not
    skimRules : list of strings
        List of skimming rules for which processors the process should listen to when deciding whether to keep an event
    logFrequency : int
        Print the event number whenever its modulus with this frequency is zero
    termLogLevel : int
        Minimum severity of log messages to print to terminal: 0 (debug) - 4 (fatal)
    fileLogLevel : int
        Minimum severity of log messages to print to file: 0 (debug) - 4 (fatal)
    logFileName : str
        File to print log messages to, won't setup file logging if this parameter is not set
    conditionsGlobalTag : str
        Global tag for the current generation of conditions
    conditionsObjectProviders : list of ConditionsObjectProviders
        List of the sources of calibration and conditions information
    randomNumberSeedService : RandomNumberSeedService
        conditions object that provides random number seeds in a deterministic way

    See Also
    --------
    Producer : one type of event processor
    Analyzer : the other type of event processor
    """

    lastProcess=None
    
    def __init__(self, passName):

        if ( Process.lastProcess is not None ) :
            raise Exception( "Process object is already created! You can only create one Process object in a script." )

        self.passName=passName
        self.maxEvents=-1
        self.maxTriesPerEvent=1
        self.run=-1
        self.inputFiles=[]
        self.outputFiles=[]
        self.sequence=[]
        self.keep=[]
        self.libraries=[]
        self.skimDefaultIsKeep=True
        self.skimRules=[]
        self.logFrequency=-1
        self.termLogLevel=2 #warnings and above
        self.fileLogLevel=0 #print all messages
        self.logFileName='' #won't setup log file
        self.compressionSetting=9
        self.histogramFile=''
        self.conditionsGlobalTag='Default'
        self.conditionsObjectProviders=[]
        self.tree_name = 'LDMX_Events'
        Process.lastProcess=self

        # needs lastProcess defined to self-register
        self.randomNumberSeedService=RandomNumberSeedService()

    def addLibrary(lib) :
        """Add a library to the list of dynamically loaded libraries

        A process object must already have been created.

        Parameters
        ----------
        lib : str
            name of library to load 

        Warnings
        --------
        - Will exit the script if a process object hasn't been defined yet.

        Examples
        --------
            addLibrary( 'libSimCore.so' )
        """

        if ( Process.lastProcess is not None ) :
            Process.lastProcess.libraries.append( lib )
        else :
            raise Exception( "No Process object defined yet! You need to create a Process before creating any EventProcessors." )
    
    def addModule(module) :
        """Add a module to the list of dynamically loaded libraries

        A process object must already have been created.

        Parameters
        ----------
        module : str
            Name of module to load as a library

        See Also
        --------
        Process.addLibrary

        Examples
        --------
        You can use this function to load a general module
            addModule('SimCore')

        With the string substitutions that are made, you can
        refer to submodules with cmake, C++, or the library
        syntax. The following calls are all equivalent.
            addModule('Ecal/Event')
            addModule('Ecal::Event')
            addModule('Ecal_Event')
        """

        actual_module_name = module.replace('/','_').replace('::','_')
        Process.addLibrary('@CMAKE_INSTALL_PREFIX@/lib/lib%s.so'%(actual_module_name))

    def declareConditionsObjectProvider(cop):
        """Declare a conditions object provider to be loaded with the process

        A process object must already have been created.

        Parameters
        ----------
        cop : ConditionsObjectProvider
            provider to load with the process

        Warnings
        --------
        - Will exit the script if a process object hasn't been defined yet.
        - Overrides an already declared COP with the passed COP if they are equal
        """

        if ( Process.lastProcess is not None ) :

            cop.setTag(Process.lastProcess.conditionsGlobalTag)

            # check if the input COP matches one already declared
            #   if it does match, override the already declared one with the passed one
            for index, already_defined_cop in enumerate(Process.lastProcess.conditionsObjectProviders) :
                if cop == already_defined_cop :
                    Process.lastProcess.conditionsObjectProviders[index] = cop
                    return

            Process.lastProcess.conditionsObjectProviders.append( cop )
        else :
            raise Exception( "No Process object defined yet! You need to create a Process before declaring any ConditionsObjectProviders." )

    def setConditionsGlobalTag(self,tag) :
        """Set the global tag for all the ConditionsObjectProviders

        Parameters
        ----------
        tag : str
            Global generation tag to pass to all COPs
        """

        self.conditionsGlobalTag=tag
        for cop in self.conditionsObjectProviders :
            cop.setTag(tag)
            
    def skimDefaultIsSave(self):
        """Configure the process to by default keep every event."""

        self.skimDefaultIsKeep=True
        
    def skimDefaultIsDrop(self):
        """Configure the process to by default drop (not save) every event."""

        self.skimDefaultIsKeep=False

    def skimConsider(self,namePat):
        """Configure the process to listen to processors matching input.

        The list of skim rules has a rather complicated form, so it
        is better to use this helper function.

        Parameters
        ----------
        namePat : str
            Pattern for the processor instanceNames to match for the Process to listen

        Example
        -------
            ecalVeto = ldmxcfg.Producer( 'ecalVeto' , 'EcalVetoProcessor' )
            # Setup of other parameters for the veto
            p.skimConsider( 'ecalVeto' )

        See Also
        --------
        skimConsiderLabelled

        """

        self.skimRules.append(namePat)
        self.skimRules.append("")

    def skimConsiderLabelled(self,namePat,labelPat):
        """Configure the process to listen to processors matching input.

        The list of skim rules has a rather complicated form, so it
        is better to use this helper function.

        Some uses of ``setStorageHint`` in processors include a "reason"
        for the given storage hint. The input label pattern is check
        on matching this "reason" for the storage hint.

        Parameters
        ----------
        namePat : str
            Pattern for the processor instanceNames to match for the Process to listen
        labelPat : str
            Pattern for the storage hint reason to match for the Process to listen

        See Also
        --------
        skimConsider

        """
        self.skimRules.append(namePat)
        self.skimRules.append(labelPat)

    def setCompression(self,algorithm,level=9):
        """set the compression settings for any output files in this process

        We combine the compression settings here in the same way that ROOT
        does. This allows the compression settings to be passed along as
        one integer rather than two without any loss of generality.

        Look at ROOT's documentation for TFile to learn more
        about the different compression algorithms and levels available
        (as well as what integers to use). There is a summary table
        below.

        Algorithm           | int | root version
        ------------------- | --- | ------------
        ROOT global default |  0  | root-6.16
        ZLIB                |  1  | root-6.16
        LZMA                |  2  | root-6.16
        Old (ROOT 5)        |  3  | root-6.16
        LZ4                 |  4  | root-6.16
        ZSTD                |  5  | root-6.20

        Level 0 : no compression is applied
        Level 9 : maximum amount of compression available from algorithm

        Parameters
        ----------
        algorithm : int
            flag for the algorithm to use
        level : int
            flag for the level of compression to use
        """

        self.compressionSetting = algorithm*100 + level

    def inputDir(self, indir) :
        """Scan the input directory and make a list of input root files to read from it

        Lists all files ending in '.root' in the input directory (not recursive).
        Extends the inputFiles list by these files.

        Parameters
        ----------
        indir : str
            Path to directory of event files to read in
        """

        import os
        fullPathDir = os.path.realpath(indir)
        self.inputFiles.extend([ os.path.join(fullPathDir,f) 
                for f in os.listdir(fullPathDir) 
                if os.path.isfile(os.path.join(fullPathDir,f)) and f.endswith('.root') 
                ])

    def parameterDump(self) :
        """Recursively extract all configuration parameters for this process

        Only includes objects somehow attached to the process.
        """

        keys_to_skip = [ 'histograms' , 'libraries' ]

        from LDMX.SimCore import simcfg
        from LDMX.Framework import histogram as h

        def extract(obj):
            """Extract the parameter from the input object"""

            if isinstance(obj,list) :
                return [ extract(o) for o in obj ]
            elif hasattr(obj,'__dict__') :
                params = dict()
                for k in obj.__dict__ :
                    if k not in keys_to_skip :
                        params[k] = extract(obj.__dict__[k])
                return params
            else :
                return obj

        return extract(self)


    def pause(self) :
        """Print this Process and wait for user confirmation to continue

        Prints the process through the print function, and then
        waits for the user to press Enter to continue.
        """

        print(self)
        input("Press Enter to continue...")

    def __str__(self):
        """Stringify this object into a human readable, helpful form.

        This function creates a very large, multi-line string that reports (almost) all of the important
        details of this configured process.

        Returns
        -------
        str
            A human-readable, multi-line description of this process object
        """

        msg = "Process with pass name '%s'"%(self.passName)
        if (self.run>0): msg += "\n using run number %d"%(self.run)
        if (self.maxEvents>0): msg += "\n Maximum events to process: %d"%(self.maxEvents)
        else: msg += "\n No limit on maximum events to process"
        if (len(self.conditionsObjectProviders)>0):
            msg += "\n conditionsObjectProviders:\n";
            for cop in self.conditionsObjectProviders:
                msg+=str(cop)
        msg += "\n Processor sequence:"
        for proc in self.sequence:
            msg += str(proc)
        if len(self.inputFiles) > 0:
            if len(self.outputFiles)==len(self.inputFiles):
                msg += "\n Files:"
                for i in range(0,len(self.inputFiles)):
                    msg += "\n  '%s' -> '%s'"%(self.inputFiles[i],self.outputFiles[i])
            else:
                msg += "\n Input files:"
                for afile in self.inputFiles:
                    msg += '\n  ' + afile
                if len(self.outputFiles) > 0:
                    msg += "\n Output file: " + self.outputFiles[0]
        elif len(self.outputFiles) > 0:
            msg += "\n Output file: " + self.outputFiles[0]
        msg += "\n Skim rules:"
        if self.skimDefaultIsKeep: msg += "\n  Default: keep the event"
        else: msg += "\n  Default: drop the event"
        for i in range(0,len(self.skimRules)-1,2):
            if self.skimRules[i+1]=="": 
                msg += "\n  Listen to hints from processors with names matching '%s'"%(self.skimRules[i])
            else:
                msg += "\n  Listen to hints with labels matching '%s' from processors with names matching '%s'"%(self.skimRules[i+1],self.skimRules[i])
        if len(self.keep) > 0:
            msg += "\n Rules for keeping previous products:"
            for arule in self.keep:
                msg += '\n  ' + arule
        if len(self.libraries) > 0:
            msg += "\n Shared libraries to load:"
            for afile in set(self.libraries):
                msg += '\n  ' + afile

        return msg

    
