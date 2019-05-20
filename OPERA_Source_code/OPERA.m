function res=OPERA(varargin)

Version='2.3';
SubVersion='2.3-beta1';
%%
%
%        _______________________________________________________________________
%       |                                                                       |
%       |   OPERA models for physchem, environmental fate and tox properties.   |
%       |                 Version 2.3 (May 2019)                                |
%       |_______________________________________________________________________|
%
%
%OPERA is a command line application developed in Matlab providing QSAR models predictions as well as
%applicability domain and accuracy assessment. All models are built on curated data from public domain.
%Molecular descriptors are calculated using PaDEL and CDK software.
%
%
%Input:
%  -s, --SDF, --MOL, --SMI  Structure file containing the molecule(s) to be
%                           predicted. IDs will be assigned if the file does not contain molecule names.
%                           Molecular descriptors will be calculated using
%                           PaDEL software. Use V2000 SDF.
%  -d, --Descriptors        pre-calculated PaDEL descriptors in a comma delimited csv file. If the first column is not
%                           "Name" as the standard PaDEL output, molecule IDs will be assinged.
%  -fp, --fingerprints      pre-calculated descriptors using CDK2.0 in a tab delimited text file.
%   -cdk, --cdk             pre-calculated fingerprints using PaDEL in a comma delimited csv file.
%  -m, --Mat, --ascii       Matlab matrix or ascii file containing PaDEL descriptors.
%  -i, --MolID              Molecule names in csv file.
%  -t, --SaltInfo           Salt IDs to improve melting point predictions. List provided in Salts.xls
%  -l, --Labels             Descriptor labels. Necessary if the descriptor file does not contain labels
%                           or contains more than the 1444 PaDEL 2D descriptors.
%
%Output:
%  -o, --Output             Output file containing the predictions, applicability domain and accuracy
%                           information. File extension could be csv or txt. The output will contain by default:
%                           Molecule ID, predicted value (pred), Applicability domain (AD), Applicability domain index
%                           (AD_index) and accuracy estimate (Conf_index).
%  -n, --Neighbors          Add 5 nearest neighbors from training set (CAS, InCHiKeys, Observed and predicted values)
%  -O, --FullOutput         Output file containing all prediction details and used descriptors in csv format.
%  -x, --Seperate           Separate output file for each endpoint.
%
%Miscellaneous:
%  -v, --Verbose            Verbose level: 0=silent (default), 1=minimum details, %  2=full details.
%  -a, --All                All endpoints to be calculated (default).
%  -c, --Clean              Remove temporary files (generated during descriptor calculation.)
%  -LogP, -BCF...           List endpoints to be calculated (case insensitive). 'BCF'/'LogBCF','BP','LogP','MP',
%                           'VP'/'LogVP','WS', 'AOH', 'BioDeg', 'RB'/'ReadyBiodeg','HL'/'LogHL','KM'/'LogKM',
%                           'KOA','Koc'/'LogKoc', 'RT', 'pKa', 'LogD', 'CERAPP'/'ER', 'CoMPARA'/'AR', 'CATMoS/AcuteTox'.
% 							Groups of Endpoints: StrP (Structural properties), PC/Physchem, EnvFate/EF, Tox (ER, AR, AcuteTox).
%  -e, --Endpoint      		List endpoints to be calculated.
%  -h, --Help               Display this help file and exit.
%  -V, --Version            Version of the application
%
%
%
%
%Developed by:
%Kamel Mansouri
%mansourikamel@gmail.com
%
%
%For more information about the models and the data:
%[1] Mansouri, K. et al. SAR and QSAR in Env. Res. (2016). https://doi.org/10.1080/1062936X.2016.1253611
%[2] Mansouri K. et al. J Cheminform (2018) https://doi.org/10.1186/s13321-018-0263-1.
%[3] The CompTox Chemistry Dashboard (https://comptox.epa.gov/dashboard)
%[4] Williams A. J. et al. J Cheminform (2017) https://doi.org/10.1186/s13321-017-0247-6
%[5] JRC QSAR Model Database https://qsardb.jrc.ec.europa.eu/qmrf/endpoint



%%
% s = SplashScreen( 'Splashscreen', 'Splash_OPERA.gif', ...
%                         'ProgressBar', 'on', ...
%                         'ProgressPosition', 5, ...
%                         'ProgressRatio', 0.4 );
% delete(s);
%%


if nargin==0
    %
    % % % Read in your GIF file. Don't forget to read in the colour map as it is
    % required for display.
    [I, map]=imread('Splash5_OPERA.gif','Frames','all');
    
    % Create a figure to hold your splashscreen
    hfig=figure;
    set(hfig,'Menubar', 'none');
    set(hfig,'name','Please wait. Loading...','numbertitle','off');
    
    % Set a timer to dynamically update the plot every 0.1 sec
    t=timer('TimerFcn', {@timerCallbackFcn, hfig, I, map},'ExecutionMode','FixedRate','Period',0.1);
    
    % Start the timer
    start(t);
    
    % Check path variable on Windows
    % if isdeployed
    %     [status, result] = system('echo %PATH%');
    %     %contains(result,'OPERA')
    %     if isempty(regexpi(result,'C:\Program Files\OPERA\application'))
    %       %[status, result] = system('IF EXIST "C:\Program Files\OPERA\application" set "PATH=%PATH%;C:\Program Files\OPERA\application"');
    %
    %       setenv('PATH', [getenv('PATH') ';C:\Program Files\OPERA\application']);
    %       system('echo %PATH%')
    %
    %     end
    % end
    
    % Do your stuff here
    for j=1:5
        pause(0.5);
    end
    %train=load ('OPERA_models.mat', '-mat');
    
    % Clean-up
    stop(t);
    delete(t);
    delete(hfig);
    %
    if ispc
        type('intro_w.txt')
    else
        type('intro.txt')
    end
    help=1;
    
else
    
    FileOut=['OPERA',Version,'_Pred.csv'];
    verbose=0;InputMatrix=0;
    importedNames=0;
    importedLabels=0;
    input=0;
    InputDescPadel=0;
    inputFP=0;
    inputCDK=0;
    structure=0;
    clean=0;
    printtDesc=0;
    sep=0;
    all=1;
    prop={};
    salt=0;
    help=0;
    e=0;
    i=1;
    neighbors=0;
    fp=0;
    cdk=0;
    InputDesc={};
    InputDescFP={};
    InputDescCDK={};
    tox=0;
    pc=0;
    ef=0;
    adme=0;
    exp=0;
    
    
    
    %if nargin>0
    while i<=length(varargin)
        if  strcmpi('--descriptors',varargin{i})|| strcmpi('-d',varargin{i})|| strcmpi('--desc',varargin{i})
            InputDesc=varargin{i+1};
            InputDescPadel=1;
            input=1;
            i=i+2;
            continue
        elseif strcmpi('--fingerprints',varargin{i})|| strcmpi('-fp',varargin{i})
            InputDescFP=varargin{i+1};
            inputFP=1;
            i=i+2;
        elseif strcmpi('--cdk',varargin{i})|| strcmpi('-cdk',varargin{i})
            InputDescCDK=varargin{i+1};
            inputCDK=1;
            i=i+2;
        elseif strcmpi('--mat',varargin{i}) || strcmpi('--ascii',varargin{i})|| strcmp('-m',varargin{i})
            InputDesc=varargin{i+1};
            InputMatrix=1;
            input=1;
            i=i+2;
            continue
        elseif strcmpi('--MolID',varargin{i})|| strcmp('-i',varargin{i})
            MolID=varargin{i+1};
            MoleculeNames=importfile(MolID);
            importedNames=1;
            i=i+2;
            continue
        elseif  strcmpi('--labels',varargin{i}) || strcmp('-l',varargin{i})
            labels=varargin{i+1};
            Xlabels=importfile(labels);
            importedLabels=1;
            i=i+2;
            continue
        elseif strcmpi('--structure',varargin{i}) || strcmpi('--sdf',varargin{i})|| strcmpi('--smiles',varargin{i})|| strcmpi('--smi',varargin{i})|| strcmpi('-s',varargin{i})|| strcmpi('--mol',varargin{i})
            StructureFile=varargin{i+1};
            %InputDesc='PadelDesc.csv';
            structure =1;
            input=1;
            %clean=1;
            i=i+2;
            continue
        elseif strcmpi('--out',varargin{i}) || strcmpi('--fullOutput',varargin{i})|| strcmpi('-o',varargin{i})
            FileOut=varargin{i+1};
            if strcmpi('--fullOutput',varargin{i})|| strcmp('-O',varargin{i})
                printtDesc=1;
                neighbors=1;
            end
            i=i+2;
            continue
        elseif strcmp('-v',varargin{i}) || strcmpi('--verbose',varargin{i})
            verbose=varargin{i+1};
            if ischar(verbose)
                verbose=str2double(verbose);
            end
            i=i+2;
            continue
        elseif strcmpi('--Clean',varargin{i})|| strcmp('-c',varargin{i})
            clean=1;
            i=i+1;
            continue
        elseif strcmpi('--Neighbors',varargin{i})|| strcmp('-n',varargin{i})
            neighbors=1;
            i=i+1;
            continue
        elseif strcmpi('--salt',varargin{i}) || strcmpi('--saltInfo',varargin{i})|| strcmp('-t',varargin{i})
            salt=1;
            FileSalt=varargin{i+1};
            i=i+2;
            continue
        elseif strcmpi('--sep',varargin{i}) || strcmpi('--separate',varargin{i})|| strcmp('-x',varargin{i})
            sep=1;
            i=i+1;
            continue
        elseif strcmpi('-All',varargin{i})|| strcmp('-a',varargin{i})
            all=1;
            fp=1;
            %InputDescFP='PadelFP.csv';
            cdk=1;
            %InputDescCDK='CDKDesc.csv';
            i=i+1;
            continue
        elseif strcmp('-e',varargin{i})|| strcmpi('--endpoint',varargin{i})
            all=0;
            e=1;
            i=i+1;
            continue
        elseif strcmpi('StrP',varargin{i}) || strcmpi('BCF',varargin{i}) || strcmpi('BP',varargin{i})|| strcmpi('LogP',varargin{i})|| strcmpi('MP',varargin{i})|| strcmpi('VP',varargin{i})|| strcmpi('WS',varargin{i})...
                || strcmpi('LogWS',varargin{i})|| strcmpi('LogVP',varargin{i})|| strcmpi('LogBCF',varargin{i})|| strcmpi('AOH',varargin{i})|| strcmpi('BioHC',varargin{i})...
                || strcmpi('Biowin',varargin{i})|| strcmpi('RB',varargin{i})|| strcmpi('HL',varargin{i})|| strcmpi('KM',varargin{i})|| strcmpi('KOA',varargin{i})...
                || strcmpi('KOC',varargin{i})|| strcmpi('LogKOC',varargin{i})|| strcmpi('LogKM',varargin{i})|| strcmpi('LogHL',varargin{i})|| strcmpi('BioDeg',varargin{i})|| strcmpi('AOH',varargin{i})...
                || strcmpi('ReadyBiodeg',varargin{i})|| strcmpi('RT',varargin{i})|| strcmpi('Rbiodeg',varargin{i})||strcmpi('BioHL',varargin{i})||strcmpi('BioDegHL',varargin{i})||strcmpi('pka',varargin{i})||strcmpi('LogD',varargin{i})||strcmpi('EnvFate',varargin{i})||strcmpi('EF',varargin{i})...
                ||strcmpi('ER',varargin{i})||strcmpi('CERAPP',varargin{i})||strcmpi('AR',varargin{i})||strcmpi('CoMPARA',varargin{i})||strcmpi('AcuteTox',varargin{i})||strcmpi('CATMoS',varargin{i})||strcmpi('Tox',varargin{i})||strcmpi('PhysChem',varargin{i})||strcmpi('PC',varargin{i})...
                ||strcmpi('FuB',varargin{i})||strcmpi('FU',varargin{i})||strcmpi('ADME',varargin{i})||strcmpi('Clint',varargin{i})||strcmpi('Cl',varargin{i})
            if strcmpi('pka',varargin{i})||strcmpi('LogD',varargin{i})||strcmpi('PhysChem',varargin{i})||strcmpi('PC',varargin{i})
                fp=1;
                
            elseif strcmpi('ER',varargin{i})||strcmpi('CERAPP',varargin{i})||strcmpi('AR',varargin{i})||strcmpi('CoMPARA',varargin{i})||strcmpi('AcuteTox',varargin{i})||strcmpi('CATMoS',varargin{i})||strcmpi('Tox',varargin{i})...
                    ||strcmpi('FuB',varargin{i})||strcmpi('FU',varargin{i})||strcmpi('ADME',varargin{i})||strcmpi('Clint',varargin{i})||strcmpi('Cl',varargin{i})
                cdk=1;
                
            end
            
            if e==1
                if strcmpi('Tox',varargin{i})
                    prop=[prop, 'CERAPP','CoMPARA', 'CATMoS'];
                    all=0;
                    tox=1;
                elseif strcmpi('PhysChem',varargin{i})||strcmpi('PC',varargin{i})
                    prop=[prop, 'BP','LogP','MP','VP','WS', 'HL', 'KOA', 'RT','pKa', 'LogD'];
                    all=0;
                    pc=1;
                elseif strcmpi('EnvFate',varargin{i})||strcmpi('EF',varargin{i})
                    prop=[prop, 'LogBCF', 'AOH', 'BioDeg', 'RBioDeg','KM','KOC'];
                    all=0;
                    ef=1;
                elseif strcmpi('ADME',varargin{i})
                    prop=[prop, 'FuB', 'Clint'];
                    all=0;
                    adme=1;
                else
                    prop=[prop varargin{i}];
                    all=0;
                end
            else
                error('Check input arguments or type -h, --help for more info.')
            end
            i=i+1;
            continue
            
        elseif strcmpi('-StrP',varargin{i}) ||strcmpi('-BCF',varargin{i}) || strcmpi('-BP',varargin{i})|| strcmpi('-LogP',varargin{i})|| strcmpi('-MP',varargin{i})|| strcmpi('-VP',varargin{i})|| strcmpi('-WS',varargin{i})...
                || strcmpi('-LogWS',varargin{i})|| strcmpi('-LogVP',varargin{i})|| strcmpi('-LogBCF',varargin{i})|| strcmpi('-AOH',varargin{i})|| strcmpi('-BioHC',varargin{i})...
                || strcmpi('-Biowin',varargin{i})|| strcmpi('-RB',varargin{i})|| strcmpi('-HL',varargin{i})|| strcmpi('-KM',varargin{i})|| strcmpi('-KOA',varargin{i})...
                || strcmpi('-KOC',varargin{i})|| strcmpi('-LogKOC',varargin{i})|| strcmpi('-LogKM',varargin{i})|| strcmpi('-LogHL',varargin{i})|| strcmpi('-BioDeg',varargin{i})|| strcmpi('-AOH',varargin{i})...
                || strcmpi('-ReadyBiodeg',varargin{i})|| strcmpi('-RT',varargin{i})|| strcmpi('-Rbiodeg',varargin{i})||strcmpi('-BioHL',varargin{i})||strcmpi('-BioDegHL',varargin{i})||strcmpi('-pka',varargin{i})||strcmpi('-LogD',varargin{i})||strcmpi('-EnvFate',varargin{i})||strcmpi('-EF',varargin{i})...
                ||strcmpi('-ER',varargin{i})||strcmpi('-CERAPP',varargin{i})||strcmpi('-AR',varargin{i})||strcmpi('-CoMPARA',varargin{i})||strcmpi('-AcuteTox',varargin{i})||strcmpi('-CATMoS',varargin{i})||strcmpi('-Tox',varargin{i})||strcmpi('-PhysChem',varargin{i})||strcmpi('-PC',varargin{i})...
                ||strcmpi('-FuB',varargin{i})||strcmpi('-FU',varargin{i})||strcmpi('-ADME',varargin{i})||strcmpi('-Clint',varargin{i})||strcmpi('-Cl',varargin{i})
            if  strcmpi('-pka',varargin{i})||strcmpi('-LogD',varargin{i})||strcmpi('-PhysChem',varargin{i})||strcmpi('-PC',varargin{i})
                fp=1;
                
            elseif strcmpi('-ER',varargin{i})||strcmpi('-CERAPP',varargin{i})||strcmpi('-AR',varargin{i})||strcmpi('-CoMPARA',varargin{i})||strcmpi('-AcuteTox',varargin{i})||strcmpi('-CATMoS',varargin{i})||strcmpi('-Tox',varargin{i})...
                    ||strcmpi('-FuB',varargin{i})||strcmpi('-FU',varargin{i})||strcmpi('-ADME',varargin{i})||strcmpi('-Clint',varargin{i})||strcmpi('-Cl',varargin{i})
                cdk=1;
                
            end
            if strcmpi('-Tox',varargin{i})
                prop=[prop, 'CERAPP', 'CoMPARA','CATMoS'];
                all=0;
                tox=1;
            elseif strcmpi('-PhysChem',varargin{i})||strcmpi('-PC',varargin{i})
                prop=[prop, 'BP','LogP','MP','VP','WS', 'HL','KOA','RT','pKa', 'LogD'];
                all=0;
                pc=1;
            elseif strcmpi('-EnvFate',varargin{i})||strcmpi('-EF',varargin{i})
                prop=[prop, 'LogBCF', 'AOH', 'BioDeg', 'RBioDeg','KM','KOC'];
                all=0;
                ef=1;
            elseif strcmpi('-ADME',varargin{i})
                prop=[prop, 'FuB', 'Clint'];
                all=0;
                adme=1;
            else
                all=0;
                prop=[prop strrep(varargin{i},'-','')];
            end
            i=i+1;
            continue
        elseif strcmpi('--help',varargin{i})|| strcmp('-h',varargin{i})
            if ispc
                type('help_w.txt')
            else
                type('help.txt')
            end
            help=1;
            i=i+1;
            continue
        elseif strcmp('-V',varargin{i})|| strcmpi('--version',varargin{i})
            fprintf(1,'Version %s.\n',SubVersion);
            help=1;
            i=i+1;
            continue
        elseif strcmp('-exp',varargin{i})
            exp=1;
            i=i+1;
            continue
        else
            error('Check input arguments or type -h, --help for more info.')
            
        end
        
    end
    
    % If no splash
    %  else
    % % error('MyComponent:incorrectType',...
    % %    'Not enough arguments. \nUsage: OPERA [OPTION]... <Input> <output>... \nType -h, --help for more info.')
    %
    %
    % type('intro_w.txt')
    % help=1;
    
    %fprintf(2,'Not enough arguments \n');
    %return
end


if help==1
    %     return
    %     %('End help file!')
    % end
    % % else
    
    res=0;
else
    
    if verbose==0 || isdeployed
        warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    else
        warning('on','MATLAB:table:ModifiedAndSavedVarnames');
    end
    
    train=load ('OPERA_models.mat', '-mat');
    if importedLabels==0
        Xlabels=train.labels;
        XlabelsFP=train.labels_fp;
    end
    
    
    if  structure ==0
        if input==0
            error('No structure file or a comma delimited file with generated PaDEL descirptors. Usage: OPERA [OPTION]... <Input> <output>... Type -h, --help for more info.')
            %   fprintf(2,'You must at least enter an input file \n');
            %   return
        end
        if fp==1 && inputFP==0
            error('No structure file or a comma delimited file with generated fingerprints. Usage: OPERA [OPTION]... <Input> <output>... Type -h, --help for more info.')
        end
        if cdk==1 && inputCDK==0
            error('No structure file or a tab delimited file with calculated CDK2.0 descriptors. Usage: OPERA [OPTION]... <Input> <output>... Type -h, --help for more info.')
        end
    end
    
    
    if all==1
        prop= {'StrP','LogBCF','BP','LogP','MP','VP','WS', 'AOH', 'BioDeg', 'RBioDeg','HL','KM','KOA','KOC','RT','pKa', 'LogD', 'CERAPP', 'FuB','Clint', 'CoMPARA', 'CATMoS'};
        if verbose >0
            fprintf(1,'\n All properties will be calculated: \nGeneral structural properties, Physchem, Env. fate, ADME and Tox Endpoints (CERAPP, CoMPARA and CATMoS)  \n');
        end
        fp=1;
        cdk=1;
        
    else
        if verbose >0
            if size(prop(:),1)>1
                endpoints=strjoin(prop(1:size(prop(:),1)-1),',  ');
                fprintf(1,'\n Endpoints to be calculated: \n %s and %s\n',upper(endpoints),upper(prop{end}));
            else
                fprintf(1,'\n Endpoint to be calculated: %s\n',upper(prop{:}));
            end
            
        end
    end
    
    
    %if structure==1 && (InputDescPadel==0||(fp==1 && inputFP==0))
    
    if ispc
        installdir=fullfile('C:','Program Files','OPERA','application');
    else
        installdir=fullfile('usr','local','bin','OPERA','application');
    end
    if ~exist(fullfile(installdir,'padel-full-1.00.jar'),'file')
        
        if isdeployed
            currentDir=ctfroot;
            
            if exist(fullfile(currentDir,'OPERA_installdir.txt'),'file')
                fid  = fopen(fullfile(currentDir,'OPERA_installdir.txt'),'r');
                installdir=strip(fread(fid,'*char')');
                fclose(fid);
                
                if ~exist(fullfile(installdir,'padel-full-1.00.jar'),'file')
                    error('Default install folder was changed during installation. Update OPERA_installdir.txt file in %s', currentDir);
                end
            else
                fid  = fopen(fullfile(currentDir,'OPERA_installdir.txt'),'w');
                fprintf(fid,'%s',installdir);
                fclose(fid);
                error('Default install folder was changed during installation. Update OPERA_installdir.txt file in %s', currentDir);
            end
            
        else
            currentDir = pwd;
            if exist(fullfile(currentDir,'OPERA_installdir.txt'),'file')
                fid  = fopen(fullfile(currentDir,'OPERA_installdir.txt'),'r');
                installdir=strip(fread(fid,'*char')');
                fclose(fid);
                %currentDir = pwd;
                if ~exist(fullfile(installdir,'padel-full-1.00.jar'),'file')
                    error(['Default install folder was changed during installation. Update OPERA_installdir.txt file in ', currentDir]);
                end
            else
                fid  = fopen(fullfile(currentDir,'OPERA_installdir.txt'),'w');
                fprintf(fid,'%s',installdir);
                fclose(fid);
                error(['Default install folder was changed during installation. Update OPERA_installdir.txt file in ', currentDir]);
            end
        end
    end
    
        %---------Output file---------
    %errmsg='Cannot write to output file \n';
    ext=FileOut(length(FileOut)-3:end);
    if sep==1
        outputname=cell(size(prop));
        %output=zeros(size(prop));
        
        FileOut=strrep(FileOut,ext,'');
        
        for i=1:length(prop)
            %         FileOut(i)=[FileOut prop(i) ext]
            outputname{i}=strrep(strjoin([FileOut '_' prop(i) ext]),' ', '');
            [output(i),errmsg]=fopen(outputname{i},'w');
        end
        FileOut=outputname;
    else
        
        [output,errmsg]=fopen(FileOut,'w');
    end
    
    
    if verbose>0 && ~isempty(errmsg)
        disp('Output file')
        error(errmsg)
        % disp(errmsg);
        %  return
    end
    %-----------------------------
    
    
    %Start input Matrix
    if InputMatrix==1
        if verbose> 0
            disp('Loading matrix of descriptors...');
        end
        load(InputDesc);
        Xin=eval(InputDesc(1:length(InputDesc)-4));
        if importedNames==0 && size(Xin,1)==1444
            %MoleculeNames=num2cell(1:1:size(Xin,1))';
            for i=1:size(Xin,1)
                MoleculeNames{i,1}=strcat('AUTOGEN_',num2str(i));
            end
            %         if verbose>0
            %             disp(' default PaDEL descriptor names considered...\n');
            %         end
        else
            for i=1:size(Xin,1)
                MoleculeNames{i,1}=strcat('AUTOGEN_',num2str(Xin(i)));
            end
            Xin(:,1)=[];
        end
     %End input Matrix
    else
        if structure==1
            if strcmpi(StructureFile(length(StructureFile)-3:end),'.smi')||strcmpi(StructureFile(length(StructureFile)-3:end),'.txt')
                if ispc
                    [~, numStruct] = system(['FINDSTR /R /N "^.*" ', strcat('"',char(StructureFile),'"'),' | FIND /C ":"']);%win
                else
                    [~, numStruct] = system(['cat ', strcat('"',char(StructureFile),'"'),' | sed "/^\s*$/d" | wc -l ']); %linux
                end
                
            elseif strcmpi(StructureFile(length(StructureFile)-3:end),'.sdf')
                if ispc
                    [~, numStruct] = system(['FINDSTR /R /N "^.*$$$$" ', strcat('"',char(StructureFile),'"') ,' | FIND /C ":"']);%win
                else
                    [~, numStruct] = system(['grep -F "$" ', strcat('"',char(StructureFile),'"'), ' | wc -l']); %linux
                end
                
            elseif strcmpi(StructureFile(length(StructureFile)-3:end),'.mol')
                
                numStruct='1';
            else
                error('Check input file');
                
            end
            
            if strcmpi(StructureFile(length(StructureFile)-3:end),'.txt')
                if verbose >0
                    fprintf(1,'Found IDs in input text file: %d.\n',str2double(numStruct));
                end
                fid = fopen(StructureFile);
                indic = 1;
                while 1
                    tline = fgetl(fid);
                    if ~ischar(tline)
                        break
                    end
                    strings{indic}=tline;
                    indic = indic + 1;
                end
                fclose(fid);
                StructureFile=strcat(StructureFile(1:length(StructureFile)-3),'smi');
                fileID = fopen(StructureFile, 'w');
                f=0;
                for i=1:length(strings)
                    [La,Lb] = ismember(strings{i},train.DSSToxQSARr{:,2});
                    if La
                        f=f+1;
                        fprintf(fileID,'%s\t%s\n',train.DSSToxQSARr{Lb,1},strings{i});
                    else
                        [La,Lb] = ismember(strings{i},train.DSSToxQSARr{:,3});
                        if La
                            f=f+1;
                            fprintf(fileID,'%s\t%s\n',train.DSSToxQSARr{Lb,1},strings{i});
                        else
                            [La,Lb] = ismember(strings{i},train.DSSToxQSARr{:,4});
                            if La
                                f=f+1;
                                fprintf(fileID,'%s\t%s\n',train.DSSToxQSARr{Lb,1},strings{i});
                            else
                                [La,Lb] = ismember(strings{i},train.DSSToxQSARr{:,5});
                                if La
                                    f=f+1;
                                    fprintf(fileID,'%s\t%s\n',train.DSSToxQSARr{Lb,1},strings{i});
                                end
                            end
                        end
                    end
                end
                fclose(fileID);
                if verbose >0
                    fprintf(1,'Found structures based on provided IDs: %d.\n',f);
                end
                numStruct=num2str(f);
            end
        end
        %========== Molecular Descriptors ==========
        
        %Calculating PaDEL descriptors...
        
        if verbose >0
            fprintf(1,'\n========== Molecular Descriptors ==========\n');
        end
        
        if structure==1 && InputDescPadel==0
            InputDesc=strcat(StructureFile(1:length(StructureFile)-4),'_PadelDesc.csv');
            if verbose >0
                fprintf(1,'Loaded structures: %d.\n',str2double(numStruct));
                
                fprintf(1,'PaDEL calculating 2D descriptors...\n');
                if verbose ==1
                    [statusDesc,cmdoutDesc] =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                        ' -2d -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir ' strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDesc,'"')...
                        ' > ' strcat('"','PaDELlogfile.log','"')]);
                else
                    statusDesc =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                        ' -2d -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir ' strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDesc,'"')]);
                end
                fprintf(1,'PaDEL descriptors calculated for: ');
                
                if ispc
                    [~, numlines] = system(['FINDSTR /R /N "^.*" ',InputDesc,' | FIND /C ":"']); %win
                else
                    [~, numlines] = system( ['wc -l ', InputDesc] ); %linux
                end
                numlines=str2double(strrep(numlines,InputDesc,''))-1;
                fprintf(1, '%d molecules.\n',numlines);
                if numlines < str2double(numStruct)
                    error('PaDEL descriptors failed. Check input structures!');
                end
                if statusDesc~=0 && ~isempty(cmdoutDesc)
                    disp(cmdoutDesc);
                end
            else
                [~,~] =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                    ' -2d -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir ' strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDesc,'"')...
                    ' > ' strcat('"','PaDELlogfile.log','"')]);
                if ispc
                    [~, numlines] = system(['FINDSTR /R /N "^.*" ',InputDesc,' | FIND /C ":"']); %win
                else
                    [~, numlines] = system( ['wc -l ', InputDesc] ); %linux
                end
                numlines=str2double(strrep(numlines,InputDesc,''))-1;
                if numlines < str2double(numStruct)
                    error('PaDEL descriptors failed. Check input structures!');
                end
            end
            
        end
        
        
        
        %if statusDesc==0 && verbose>0 && isempty(cmdoutDesc)
        %fprintf(1,'--------------------------------------------------------\n');
        
        %Windows OS: store the below two lines in countlines.pl
        %while (<>) {};
        %print $.,"\n";
        %Then to make a matlab call to count the lines for file XYZ.csv
        %numlines = str2num( perl('countlines.pl', 'XYZ.csv') );

        %end
        
        if verbose> 0
            disp('Loading of PaDEL descriptors file...');
        end
        
        %Xin=dataset('File',InputDesc,'delimiter',',');
        try
            Xin=readtable(InputDesc,'delimiter',',');
        catch ME
            if strcmp(ME.identifier,'MATLAB:readtable:OpenFailed')
                error('Unable to open descriptors file');
            else
                error(ME.message);
                return;
            end
        end
        %Xlabels=Xin.Properties.VarNames;
        Xlabels=Xin.Properties.VariableNames;
        %Xin=dataset2table(Xin);
        
        
        if size(Xin,1)==0 || size(Xin,2)==0
            error('Check input file and re-run to calculate descriptors')
            %     elseif verbose>0
            %         fprintf(1,'The number of input molecules is: %d \n',size(Xin,1));
            
        end
        
        
        if strcmpi(Xlabels{1},'Name') || strcmpi(Xlabels{1},'MoleculeID') || strcmpi(Xlabels{1},'Molecule')
            
            if verbose> 1
                disp('Molecule names found in input file(s).');
            end
            Xin.Properties.VariableNames{1}='Name';
            Xlabels=Xlabels(2:end);
            Names=Xin.Name;
            if isnumeric(Names) && strcmpi(ext,'.txt')
                
                for i=1:size(Xin,1)
                    MoleculeNames{i,1}=strcat('AUTOGEN_',num2str(Names(i)));
                end
            else
                MoleculeNames=Names;
                
                %MoleculeID=num2cell(MoleculeNames);
                
            end
            %Xin=Xin{:,2:end};
            Xin=Xin(:,2:end);
        else
            
            if verbose> 1
                disp('Molecule names not found in input file. Generated IDs will be assigned.');
            end
            %Xin=Xin{:,:};
            %Xin=Xin(:,:);
            
            for i=1:size(Xin,1)
                MoleculeNames{i,1}=strcat('AUTOGEN_',num2str(i));
            end
        end
        if size(Xin,1)==0 || size(Xin,2)==0
            error('Empty descriptors file!');
        end
        i=1;
        Temp=zeros(size(Xin));
        if verbose> 0
            disp('Checking loaded variables.');
        end
        while i<=length(Xlabels)
            if cellfun(@ischar,table2cell(Xin(1,i)))
                Temp(:,i)=str2double(table2cell(Xin(:,i)));
            else
                
                Temp(:,i)=Xin{:,i};
            end
            i=i+1;
        end
        if verbose> 0
            disp(['Loaded ', num2str(length(Xlabels)),' PaDEL descriptors for ', num2str(size(Xin,1)),' molecules.']);
        end
        
        clear('Xin');
        Xin=Temp;
        clear('Temp');
 
        if cdk==1
            if structure==1 && inputCDK==0
                %Bond_HA_r=Xin(:,466)./Xin(:,9);
                Amb_str=intersect(find((Xin(:,466)./Xin(:,9))>1.3),find(Xin(:,9)>50));
                if ~isempty(Amb_str)||~isempty(find(Xin(:,9)>150, 1))
                    Amb_str=num2str(Amb_str);
                    Amb_str=strjoin(num2cell(Amb_str(1:length(Amb_str))),',  ');
                    error('CDK descriptors cannot be calculated for structure(s) number: %s.',Amb_str);
                end
                
                InputDescCDK=strcat(StructureFile(1:length(StructureFile)-4),'_CDKDesc.csv');
                fprintf(1,'CDK 2.0 calculating 2D descriptors...\n');
                [statusDesc,cmdoutDesc] =system (['java -jar ' strcat('"',fullfile(installdir,'CDKDescUI-2.0.jar'),'"') ' -b -t all -o ' strcat('"',InputDescCDK,'"')...
                    ' ' strcat('"',char(StructureFile),'"') ' > ' strcat('"','CDKlogfile.log','"') ' 2> ' strcat('"','CDKerr.log','"')]);
                fprintf(1,'CDK descriptors calculated for: ');
                
                if ispc
                    [~, numlines] = system(['FINDSTR /R /N "^.*" ',InputDescCDK,' | FIND /C ":"']); %win
                else
                    [~, numlines] = system( ['wc -l ', InputDescCDK] ); %linux
                end
                numlines=str2double(strrep(numlines,InputDescCDK,''))-1;
                if verbose>0
                    fprintf(1, '%d molecules.\n',numlines);
                    if statusDesc~=0 && ~isempty(cmdoutDesc)
                        disp(cmdoutDesc);
                    end
                end
                if numlines < str2double(numStruct)
                    error('CDK descriptors failed. Check input structures!');
                end
                
                
            end
            if verbose> 0
                disp('Loading of CDK descriptors file...');
            end
            try
                XinCDK=readtable(InputDescCDK,'delimiter','\t');
            catch ME
                if strcmp(ME.identifier,'MATLAB:readtable:OpenFailed')
                    error('Unable to open descriptors file');
                else
                    error(ME.message);
                    return;
                end
            end
            if size(XinCDK,1)==0 || size(XinCDK,2)==0
                error('Empty descriptors file!');
            end
            XinCDK(:,1)=[];
            if size(XinCDK,1)~=size(Xin,1)
                error('Mismatch between PaDEL and CDK descriptors files')
            elseif strcmpi(XinCDK.Properties.VariableNames(end),'Zagreb')
                XlabelsCDK=XinCDK.Properties.VariableNames;
            elseif strcmpi(XinCDK.Properties.VariableNames(end),'nAcid')
                XinCDK=XinCDK(:,train.reorder_CDK);
                XlabelsCDK=XinCDK.Properties.VariableNames;
            else
                error('Check or recalculate CDK descriptors');
            end
            if size(XinCDK,1)==size(Xin,1)
                %fprintf(1,'The number of input molecules is: %d \n',size(XinCDK,1));
                
                i=1;
                Temp=zeros(size(XinCDK));
                if verbose> 0
                    disp('Checking loaded variables.');
                end
                while i<=length(XlabelsCDK)
                    if cellfun(@ischar,table2cell(XinCDK(1,i)))
                        Temp(:,i)=str2double(table2cell(XinCDK(:,i)));
                    else
                        Temp(:,i)=XinCDK{:,i};
                    end
                    i=i+1;
                end
                if verbose> 0
                    %disp(['The number of loaded CDK descriptors is: ', num2str(length(XlabelsCDK))]);
                    disp(['Loaded ', num2str(length(XlabelsCDK)),' PaDEL descriptors for ', num2str(size(XinCDK,1)),' molecules.']);
                    
                end
                clear('XinCDK');
                XinCDK=Temp;
                clear('Temp');
            end
            
        end
        if fp==1
            if structure==1 && inputFP==0
                InputDescFP=strcat(StructureFile(1:length(StructureFile)-4),'_PadelFP.csv');
                if verbose >0
                    fprintf(1,'PaDEL generating fingerprints...\n');
                    if verbose ==1
                        [statusDesc,cmdoutDesc] =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                            ' -fingerprints -descriptortypes ' strcat('"',fullfile(installdir,'desc_fp.xml'),'"') ' -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir '...
                            strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDescFP,'"') ' > ' strcat('"','PaDELlogfileFP.log','"')]);
                    else
                        [statusDesc,cmdoutDesc] =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                            ' -fingerprints -descriptortypes ' strcat('"',fullfile(installdir,'desc_fp.xml'),'"') ' -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir '...
                            strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDescFP,'"')]);
                    end
                    fprintf(1,'PaDEL fingerprints generated for: ');
                    
                    if ispc
                        [~, numlines] = system(['FINDSTR /R /N "^.*" ',InputDescFP,' | FIND /C ":"']); %win
                    else
                        [~, numlines] = system( ['wc -l ', InputDescFP] ); %linux
                    end
                    numlines=str2double(strrep(numlines,InputDescFP,''))-1;
                    fprintf(1, '%d molecules.\n',numlines);
                    if numlines < str2double(numStruct)
                        error('PaDEL fingerprints failed. Check input structures!');
                    end
                    if statusDesc~=0 && ~isempty(cmdoutDesc)
                        disp(cmdoutDesc);
                    end
                    
                else
                    [~,~] =system (['java -Djava.awt.headless=true -jar ' strcat('"',fullfile(installdir,'padel-full-1.00.jar'),'"')...
                        ' -fingerprints -descriptortypes ' strcat('"',fullfile(installdir,'desc_fp.xml'),'"') ' -removesalt -standardizenitro -detectaromaticity -retainorder -maxruntime 60000 -dir '...
                        strcat('"',char(StructureFile),'"') ' -file ' strcat('"',InputDescFP,'"') ' > ' strcat('"','PaDELlogfileFP.log','"')]);
                    if ispc
                        [~, numlines] = system('FINDSTR /R /N "^.*" PadelFP.csv | FIND /C ":"'); %win
                    else
                        [~, numlines] = system( ['wc -l ', 'PadelFP.csv'] ); %linux
                    end
                    numlines=str2double(strrep(numlines,' PadelFP.csv',''))-1;
                    if numlines < str2double(numStruct)
                        error('PaDEL fingerprints failed. Check input structures!');
                    end
                end
                
            end
            
            if verbose> 0
                disp('Loading of fingerprints file...');
            end
            try
                XinFP=readtable(InputDescFP,'delimiter',',');
            catch ME
                if strcmp(ME.identifier,'MATLAB:readtable:OpenFailed')
                    error('Unable to open descriptors file');
                else
                    error(ME.message);
                    return;
                end
            end
            XlabelsFP=XinFP.Properties.VariableNames;
            if size(XinFP,1)==0 || size(XinFP,2)==0
                error('Empty descriptors file!');
            end
            XinFP=XinFP(:,2:end);
            %XlabelsFP=XlabelsFP(2:end);
            if size(XinFP,1)~=size(Xin,1)
                error('Mismatch between descriptors and fingerprint files')
            elseif verbose>0
                %fprintf(1,'The number of input molecules is: %d \n',size(XinFP,1));
                %disp(['The number of loaded fingerprints bits is: ', num2str(length(XlabelsFP)-1)]);
                disp(['Loaded ', num2str(length(XlabelsFP)-1),' PaDEL fingerprints for ', num2str(size(XinFP,1)),' molecules.']);
            end
            
        end
        
        %Start SaltInfo
        if salt==1
            if verbose> 0
                disp('Reading file with salt information.');
            end
            SaltIndex=readtable(FileSalt,'delimiter',',');
            %     if strcmpi(SaltIndex{1},'Name')
            %         SaltIndex=SaltIndex(
            if size(SaltIndex,1)==size(Xin,1)
                if cellfun(@ischar,table2cell(SaltIndex(1,end)))
                    Temp(:,:)=str2double(table2cell(SaltIndex(:,end)));
                else
                    
                    Temp(:,:)=SaltIndex{:,end};
                end
                
                clear('SaltIndex');
                SaltIndex=Temp;
                clear('Temp');
                
                if verbose> 0
                    disp(['The number of molecules with salt information:', num2str(length(find(SaltIndex)))]);
                end
            else
                error('The number of compounds must be the same in both files.')
                %fprintf(2,'Number of compounds must be the same in both files. \n');
                %return
            end
            
            %res.SaltID=SaltIndex;
        end
        %End SaltInfo
        
        if clean==1 && structure==1
            delete(InputDesc);
            %delete(strcat('PadelDesc_',StructureFile(1:length(StructureFile)-3),'.csv'));
            if verbose <2
                delete('PaDELlogfile.log');
            end
            if fp==1
                delete(InputDescFP);
                %delete('PadelFP.csv');
                if verbose <2
                    delete('PaDELlogfileFP.log');
                end
            end
            if cdk==1
                delete(InputDescCDK);
                %delete('CDKDesc.csv');
                if verbose <2
                    delete('CDKlogfile.log');
                    delete('CDKerr.log');
                end
            end
            
        end
        
    end
    
    
    
    if verbose> 0 && size(prop(:),1)>=1
        fprintf(1,'\n========== Running The Models ==========\n');
    end
    
    
    % General Structural properties:
    [Lia,Locb] = ismember('strp',lower(prop));
    if Lia
        if verbose> 0 && size(prop(:),1)>1
            fprintf(1,'Generating the general structural properties...\n');
        end
        %             Desc={'MW','nAtom','nHeavyAtom','nC','nO','nN','naAromAtom','nRing','nHeteroRing','HybRatio','nRotB','nHBAcc','nHBDon','LipinskiFailures','TopoPSA','AMR','MLFER_S'};
        %              Xtest=zeros(size(Xin,1),length(Desc));
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.StrP.Desc_i);
        
        Desc={'MolWeight','nbAtoms','nbHeavyAtoms','nbC','nbO','nbN','nbAromAtom','nbRing','nbHeteroRing','Sp3Sp2HybRatio','nbRotBd','nbHBdAcc','ndHBdDon','nbLipinskiFailures','TopoPolSurfAir','MolarRefract','CombDipolPolariz'};
        Xtest=array2table(Xtest,'VariableNames',Desc);
        T=array2table(MoleculeNames,'VariableNames',{'MoleculeID'});
        T=[T Xtest];
        if sep==1
            if strcmpi(ext,'.csv')
                %T=struct2table(res);
                %                     res.Descriptors=Xtest;
                
                writetable(T,FileOut{Locb},'Delimiter',',');%,'QuoteStrings',true);
                fclose(output(Locb));
            elseif strcmpi(ext,'.txt')
                fprintf(output(Locb(1)),'\n\n\t\t\t\t\t General Structural properties... \n\n			============================================================== \n\n');
                for i=1:size(Xtest(:,1))
                    fprintf(output(Locb(1)),'\t Molecule %s:\n', MoleculeNames{i});
                    for j=1:length(Desc)
                        fprintf(output(Locb(1)),'%s= %.3f\t;\t', Desc{j},Xtest{i,j});
                    end
                    fprintf(output(Locb(1)),'\n');
                end
                
            end
            
        elseif strcmpi(ext,'.txt')
            fprintf(output,'\n\n\t\t\t\t\t General Structural properties... \n\n			============================================================== \n\n');
            for i=1:size(Xtest(:,1))
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                for j=1:length(Desc)
                    fprintf(output,'%s= %.3f\t;\t', Desc{j},Xtest{i,j});
                end
                fprintf(output,'\n');
            end
        end
        res=table2struct(T,'ToScalar',true);
        
        if sep==1
            resf.StrP=table2struct(T,'ToScalar',true);
            clear('res');
        end
        clear('T');
    end
    
    
    DescMat=[];
    DescNames={};
    
    
    %for j=1:length(prop)
    %switch lower(prop{j})
    
    if verbose> 0 && (pc||all)
        fprintf(1,'---------- PhysChem properties ----------\n');
    end
    
    %Predict LogP values
    %case {'logp'}
    [Lia,Locb] =ismember({'logp','logd'},lower(prop));
    if find(Lia)
        
        Desc=train.LogP.Desc;
        
        if verbose>0
            disp('Predicting LogP values (Log10)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0 && Lia(1)
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogP values... \n\n			==============================================================n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.LogP.Desc_i);
        
        pred = nnrpred(Xtest,train.LogP.model.set.train,train.LogP.model.set.y,train.LogP.model.set.K,train.LogP.model.set.dist_type,train.LogP.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogP_exp=NaN(size(Xtest,1),1);
        end
        res.LogP_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.LogP.model.set.train,Xtest,'auto');
        res.AD_LogP=abs(AD.inorout-1)';
        res.AD_LogP(round(pred.D,3)==0)=1;
        
        %             res.AD_index1=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index1)
        %                 res.AD_index1(isnan(res.AD_index1))=0;
        %             end
        
        %res.AD_index1=1./(1+nanmedian(pred.dc,2));
        
        
        res.AD_index_LogP=zeros(size(Xtest,1),1);
        %             res.Conf_index1=zeros(size(Xtest,1),1);
        res.Conf_index_LogP=zeros(size(Xtest,1),1);
        
        %             res.dc=pred.dc;
        %             res.w=pred.w;
        LogP_CAS_neighbor=cell(size(Xtest,1),5);
        LogP_InChiKey_neighbor=cell(size(Xtest,1),5);
        LogP_DTXSID_neighbor=cell(size(Xtest,1),5);
        LogP_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogP_Exp_neighbor=zeros(size(Xtest,1),5);
        LogP_pred_neighbor=zeros(size(Xtest,1),5);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.LogP.CAS);
                    if Li
                        res.LogP_exp(i,1)=train.LogP.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.LogP.DTXSID);
                        if Li
                            res.LogP_exp(i)=train.LogP.model.set.y(Lo);
                        end
                    end
                end
            end
            
            LogP_CAS_neighbor(i,:)=train.LogP.CAS(pred.neighbors(i,:));
            LogP_InChiKey_neighbor(i,:)=train.LogP.InChiKey(pred.neighbors(i,:));
            LogP_DTXSID_neighbor(i,:)=train.LogP.DTXSID(pred.neighbors(i,:));
            LogP_DSSTOXMPID_neighbor(i,:)=train.LogP.DSSTOXMPID(pred.neighbors(i,:));
            LogP_Exp_neighbor(i,:)=train.LogP.model.set.y(pred.neighbors(i,:));
            LogP_pred_neighbor(i,:)=train.LogP.model.yc_weighted(pred.neighbors(i,:));
            
            res.AD_index_LogP(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_LogP(i,1)=((1/(1+sqrt(((LogP_Exp_neighbor(i,:)-LogP_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_LogP(i,1))/2;
            
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogP_pred(i,1)=NaN;
                res.AD_LogP(i)=0;
                res.AD_index_LogP(i)=0;
                res.Conf_index_LogP(i,1)=0;
            end
            %res.Conf_index_LogP(i,1)=((1/(1+sqrt(((LogP_Exp_neighbor(i,:)-LogP_pred_neighbor(i,:)).^2)*pred.w(i,:)'))));
            
            
            %                 rmse=calc_reg_param(res.LogP_Exp_neighbor(i,:),res.LogP_pred_neighbor(i,:));
            %                 res.Conf_index1(i,1)=1/(1+rmse.RMSEC);
            
            %res.Conf_index(i,1)=1/(1+sqrt(sum(diag((res.LogP_Exp_neighbor(i,:)-res.LogP_pred_neighbor(i,:))*pred.w(i,:)').^2)));
            
            if neighbors==1
                res.LogP_CAS_neighbor(i,:)=LogP_CAS_neighbor(i,:);
                res.LogP_InChiKey_neighbor(i,:)=LogP_InChiKey_neighbor(i,:);
                res.LogP_DTXSID_neighbor(i,:)=LogP_DTXSID_neighbor(i,:);
                res.LogP_DSSTOXMPID_neighbor(i,:)=LogP_DSSTOXMPID_neighbor(i,:);
                res.LogP_Exp_neighbor(i,:)=LogP_Exp_neighbor(i,:);
                res.LogP_pred_neighbor(i,:)=LogP_pred_neighbor(i,:);
            end
            
            
            
            if strcmpi(ext,'.txt') && sep==1 && Lia(1)
                %res.Xtest=Xtest;
                fprintf(output(Locb(1)),'\t Molecule %s:\n', res.MoleculeID{i});
                if exp
                    fprintf(output(Locb(1)),'LogP experimental= %.3f\n', res.LogP_exp(i));
                end
                fprintf(output(Locb(1)),'LogP predicted= %.3f\n', res.LogP_pred(i));
                if res.AD_LogP(i)==1
                    fprintf(output(Locb(1)),'AD: inside\n');
                else
                    fprintf(output(Locb(1)),'AD: outside\n');
                end
                fprintf(output(Locb(1)),'AD_index= %.2f\n', res.AD_index_LogP(i));
                fprintf(output(Locb(1)),'Conf_index= %.2f\n', res.Conf_index_LogP(i));
                %CAS=strjoin(res.LogP_CAS_neighbor(i,1:5),',\t');
                %calc=strjoin(num2cell(res.LogP_pred_neighbor(i,1:5)),', ');
                %exp=strjoin(num2cell(res.LogP_Exp_neighbor(i,1:5)),', ');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(1)),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.LogP.model.set.K,res.LogP_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(1)),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.LogP.model.set.K, res.LogP_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(1)),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.LogP.model.set.K, res.LogP_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0 && Lia(1)
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n',res.MoleculeID{i});
                if exp
                    fprintf(output,'LogP experimental= %.3f\n', res.LogP_exp(i));
                end
                fprintf(output,'LogP predicted= %.3f\n', res.LogP_pred(i));
                if res.AD_LogP(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_LogP(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_LogP(i));
                %CAS=strjoin(res.LogP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.LogP.model.set.K, res.LogP_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.LogP.model.set.K, res.LogP_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.LogP.model.set.K, res.LogP_pred_neighbor(i,1:5));
                end
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv') && Lia(1)
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(1)},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(1)));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv') && Lia(1)
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
            
            
        end
        
        if sep==1
            resf.LogP=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    %Predict MP values
    [Lia,Locb] =ismember('mp',lower(prop));
    if find(Lia)
        %case 'mp'
        
        
        
        Desc=train.MP.Desc;
        
        
        if verbose>0
            disp('Predicting MP values (Deg. C)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            if salt==1
                disp('Salt information will be considered in the predictions');
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting MP values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.MP.Desc_i);
        
        AD=classical_leverage(train.MP.model.set.train,Xtest,'auto');
        
        if salt ==1
            Xtest=[Xtest SaltIndex];
            Desc=[Desc,'SaltIndex'];
            pred = nnrpred(Xtest,train.MP.model_s.set.train,train.MP.model_s.set.y,train.MP.model_s.set.K,train.MP.model_s.set.dist_type,train.MP.model_s.set.param.pret_type);
            pred.D=diag(pred.D);
            %AD=classical_leverage(train.MP.model_s.set.train,Xtest,'auto');
        else
            pred = nnrpred(Xtest,train.MP.model.set.train,train.MP.model.set.y,train.MP.model.set.K,train.MP.model.set.dist_type,train.MP.model.set.param.pret_type);
            pred.D=diag(pred.D);
            %AD=classical_leverage(train.MP.model.set.train,Xtest,'auto');
        end
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.MP_exp=NaN(size(Xtest,1),1);
        end
        res.MP_pred(:,1)=pred.y_pred_weighted;
        %AD=classical_leverage(train.MP.model.set.train,Xtest,'auto');
        res.AD_MP=abs(AD.inorout-1)';
        res.AD_MP(round(pred.D,3)==0)=1;
        
        %            res.AD_index1=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        res.AD_index_MP=zeros(size(Xtest,1),1);
        res.Conf_index_MP=zeros(size(Xtest,1),1);
        
        MP_CAS_neighbor=cell(size(Xtest,1),5);
        MP_InChiKey_neighbor=cell(size(Xtest,1),5);
        MP_DTXSID_neighbor=cell(size(Xtest,1),5);
        MP_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        MP_Exp_neighbor=zeros(size(Xtest,1),5);
        MP_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.MP.CAS);
                    if Li
                        res.MP_exp(i,1)=train.MP.model_s.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.MP.DTXSID);
                        if Li
                            res.MP_exp(i)=train.MP.model_s.set.y(Lo);
                        end
                    end
                end
            end
            MP_CAS_neighbor(i,:)=train.MP.CAS(pred.neighbors(i,:));
            MP_InChiKey_neighbor(i,:)=train.MP.InChiKey(pred.neighbors(i,:));
            MP_DTXSID_neighbor(i,:)=train.MP.DTXSID(pred.neighbors(i,:));
            MP_DSSTOXMPID_neighbor(i,:)=train.MP.DSSTOXMPID(pred.neighbors(i,:));
            if salt ==1
                MP_Exp_neighbor(i,:)=train.MP.model_s.set.y(pred.neighbors(i,:));
                MP_pred_neighbor(i,:)=train.MP.model_s.yc_weighted(pred.neighbors(i,:));
            else
                MP_Exp_neighbor(i,:)=train.MP.model.set.y(pred.neighbors(i,:));
                MP_pred_neighbor(i,:)=train.MP.model.yc_weighted(pred.neighbors(i,:));
            end
            
            %                 rmse=calc_reg_param(res.MP_Exp_neighbor(i,:),res.MP_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC/50);
            
            res.AD_index_MP(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_MP(i,1)=((1/(1+sqrt(((MP_Exp_neighbor(i,:)-MP_pred_neighbor(i,:)).^2)*pred.w(i,:)')/50))+res.AD_index_MP(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.MP_pred(i,1)=NaN;
                res.AD_MP(i)=0;
                res.AD_index_MP(i)=0;
                res.Conf_index_MP(i,1)=0;
            end
            if neighbors==1
                res.MP_CAS_neighbor(i,:)=MP_CAS_neighbor(i,:);
                res.MP_InChiKey_neighbor(i,:)=MP_InChiKey_neighbor(i,:);
                res.MP_DTXSID_neighbor(i,:)=MP_DTXSID_neighbor(i,:);
                res.MP_DSSTOXMPID_neighbor(i,:)=MP_DSSTOXMPID_neighbor(i,:);
                res.MP_Exp_neighbor(i,:)=MP_Exp_neighbor(i,:);
                res.MP_pred_neighbor(i,:)=MP_pred_neighbor(i,:);
            end
            
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb),'MP experimental= %.3f\n', res.MP_exp(i));
                end
                fprintf(output(Locb),'MP predicted= %.3f\n', res.MP_pred(i));
                if res.AD_MP(i)==1
                    fprintf(output(Locb),'AD: inside\n');
                else
                    fprintf(output(Locb),'AD: outside\n');
                end
                fprintf(output(Locb),'AD_index= %.2f\n', res.AD_index_MP(i));
                fprintf(output(Locb),'Conf_index= %.2f\n', res.Conf_index_MP(i));
                %CAS=strjoin(res.MP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.MP.model.set.K, res.MP_CAS_neighbor{i,1:5});
                    fprintf(output(Locb),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.MP.model.set.K, res.MP_Exp_neighbor(i,1:5));
                    fprintf(output(Locb),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.MP.model.set.K, res.MP_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'MP experimental= %.3f\n', res.MP_exp(i));
                end
                fprintf(output,'MP predicted= %.3f\n', res.MP_pred(i));
                if res.AD_MP(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_MP(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_MP(i));
                %CAS=strjoin(res.MP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.MP.model.set.K, res.MP_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.MP.model.set.K, res.MP_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.MP.model.set.K, res.MP_pred_neighbor(i,1:5));
                end
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.MP=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    %Predict BP values
    [Lia,Locb] =ismember('bp',lower(prop));
    if find(Lia)
        %case 'bp'
        
        
        
        Desc=train.BP.Desc;
        
        
        if verbose>0
            disp('Predicting BP values (Deg. C)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting BP values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.BP.Desc_i);
        
        
        pred = nnrpred(Xtest,train.BP.model.set.train,train.BP.model.set.y,train.BP.model.set.K,train.BP.model.set.dist_type,train.BP.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.BP_exp=NaN(size(Xtest,1),1);
        end
        res.BP_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.BP.model.set.train,Xtest,'auto');
        res.AD_BP=abs(AD.inorout-1)';
        res.AD_BP(round(pred.D,3)==0)=1;
        
        
        %             res.AD_index_BP=1./(1+nanmedian(pred.dc,2));
        
        %            res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_BP=zeros(size(Xtest,1),1);
        res.Conf_index_BP=zeros(size(Xtest,1),1);
        
        BP_CAS_neighbor=cell(size(Xtest,1),5);
        BP_InChiKey_neighbor=cell(size(Xtest,1),5);
        BP_DTXSID_neighbor=cell(size(Xtest,1),5);
        BP_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        BP_Exp_neighbor=zeros(size(Xtest,1),5);
        BP_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.BP.CAS);
                    if Li
                        res.BP_exp(i,1)=train.BP.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.BP.DTXSID);
                        if Li
                            res.BP_exp(i)=train.BP.model.set.y(Lo);
                        end
                    end
                end
            end
            BP_CAS_neighbor(i,:)=train.BP.CAS(pred.neighbors(i,:));
            BP_InChiKey_neighbor(i,:)=train.BP.InChiKey(pred.neighbors(i,:));
            BP_DTXSID_neighbor(i,:)=train.BP.DTXSID(pred.neighbors(i,:));
            BP_DSSTOXMPID_neighbor(i,:)=train.BP.DSSTOXMPID(pred.neighbors(i,:));
            BP_Exp_neighbor(i,:)=train.BP.model.set.y(pred.neighbors(i,:));
            BP_pred_neighbor(i,:)=train.BP.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(BP_Exp_neighbor(i,:),BP_pred_neighbor(i,:));
            %                 res.Conf_index_BP(i,1)=1/(1+rmse.RMSEC/50);
            
            res.AD_index_BP(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_BP(i,1)=((1/(1+sqrt(((BP_Exp_neighbor(i,:)-BP_pred_neighbor(i,:)).^2)*pred.w(i,:)')/50))+res.AD_index_BP(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.BP_pred(i,1)=NaN;
                res.AD_BP(i)=0;
                res.AD_index_BP(i)=0;
                res.Conf_index_BP(i,1)=0;
            end
            if neighbors==1
                res.BP_CAS_neighbor(i,:)=BP_CAS_neighbor(i,:);
                res.BP_InChiKey_neighbor(i,:)=BP_InChiKey_neighbor(i,:);
                res.BP_DTXSID_neighbor(i,:)=BP_DTXSID_neighbor(i,:);
                res.BP_DSSTOXMPID_neighbor(i,:)=BP_DSSTOXMPID_neighbor(i,:);
                res.BP_Exp_neighbor(i,:)=BP_Exp_neighbor(i,:);
                res.BP_pred_neighbor(i,:)=BP_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb),'BP experimental= %.3f\n', res.BP_exp(i));
                end
                fprintf(output(Locb),'BP predicted= %.3f\n', res.BP_pred(i));
                if res.AD_BP(i)==1
                    fprintf(output(Locb),'AD: inside\n');
                else
                    fprintf(output(Locb),'AD: outside\n');
                end
                fprintf(output(Locb),'AD_index= %.2f\n', res.AD_index_BP(i));
                fprintf(output(Locb),'Conf_index= %.2f\n', res.Conf_index_BP(i));
                %CAS=strjoin(res.BP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.BP.model.set.K, res.BP_CAS_neighbor{i,1:5});
                    fprintf(output(Locb),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.BP.model.set.K, res.BP_Exp_neighbor(i,1:5));
                    fprintf(output(Locb),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.BP.model.set.K, res.BP_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'BP experimental= %.3f\n', res.BP_exp(i));
                end
                fprintf(output,'BP predicted= %.3f\n', res.BP_pred(i));
                if res.AD_BP(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_BP(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_BP(i));
                %CAS=strjoin(res.BP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.BP.model.set.K, res.BP_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.BP.model.set.K, res.BP_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.BP.model.set.K, res.BP_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.BP=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    %Predict VP values
    %case {'vp' ,'logvp'}
    [Lia,Locb] =ismember({'vp','logvp'},lower(prop));
    if find(Lia)
        
        
        
        Desc=train.VP.Desc;
        
        
        if verbose>0
            disp('Predicting LogVP values (Log10 mmHg)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting VP values in Log mmHg... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.VP.Desc_i);
        
        pred = nnrpred(Xtest,train.VP.model.set.train,train.VP.model.set.y,train.VP.model.set.K,train.VP.model.set.dist_type,train.VP.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogVP_exp=NaN(size(Xtest,1),1);
        end
        res.LogVP_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.VP.model.set.train,Xtest,'auto');
        res.AD_VP=abs(AD.inorout-1)';
        res.AD_VP(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_VP=zeros(size(Xtest,1),1);
        res.Conf_index_VP=zeros(size(Xtest,1),1);
        
        LogVP_CAS_neighbor=cell(size(Xtest,1),5);
        LogVP_InChiKey_neighbor=cell(size(Xtest,1),5);
        LogVP_DTXSID_neighbor=cell(size(Xtest,1),5);
        LogVP_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogVP_Exp_neighbor=zeros(size(Xtest,1),5);
        LogVP_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.VP.CAS);
                    if Li
                        res.LogVP_exp(i,1)=train.VP.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.VP.DTXSID);
                        if Li
                            res.LogVP_exp(i)=train.VP.model.set.y(Lo);
                        end
                    end
                end
            end
            LogVP_CAS_neighbor(i,:)=train.VP.CAS(pred.neighbors(i,:));
            LogVP_InChiKey_neighbor(i,:)=train.VP.InChiKey(pred.neighbors(i,:));
            LogVP_DTXSID_neighbor(i,:)=train.VP.DTXSID(pred.neighbors(i,:));
            LogVP_DSSTOXMPID_neighbor(i,:)=train.VP.DSSTOXMPID(pred.neighbors(i,:));
            LogVP_Exp_neighbor(i,:)=train.VP.model.set.y(pred.neighbors(i,:));
            LogVP_pred_neighbor(i,:)=train.VP.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogVP_Exp_neighbor(i,:),res.LogVP_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_VP(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            res.Conf_index_VP(i,1)=((1/(1+sqrt(((LogVP_Exp_neighbor(i,:)-LogVP_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+ res.AD_index_VP(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogVP_pred(i,1)=NaN;
                res.AD_VP(i)=0;
                res.AD_index_VP(i)=0;
                res.Conf_index_VP(i,1)=0;
            end
            if neighbors==1
                res.LogVP_CAS_neighbor(i,:)=LogVP_CAS_neighbor(i,:);
                res.LogVP_InChiKey_neighbor(i,:)=LogVP_InChiKey_neighbor(i,:);
                res.LogVP_DTXSID_neighbor(i,:)=LogVP_DTXSID_neighbor(i,:);
                res.LogVP_DSSTOXMPID_neighbor(i,:)=LogVP_DSSTOXMPID_neighbor(i,:);
                res.LogVP_Exp_neighbor(i,:)=LogVP_Exp_neighbor(i,:);
                res.LogVP_pred_neighbor(i,:)=LogVP_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogVP experimental= %.3f\n', res.LogVP_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogVP predicted= %.3f\n', res.LogVP_pred(i));
                if res.AD_VP(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_VP(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_VP(i));
                %CAS=strjoin(res.LogVP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.VP.model.set.K, res.LogVP_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.VP.model.set.K, res.LogVP_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.VP.model.set.K, res.LogVP_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0

                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogVP experimental= %.3f\n', res.LogVP_exp(i));
                end
                fprintf(output,'LogVP predicted= %.3f\n', res.LogVP_pred(i));
                if res.AD_VP(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_VP(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_VP(i));
                %CAS=strjoin(res.LogVP_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.VP.model.set.K, res.LogVP_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.VP.model.set.K, res.LogVP_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.VP.model.set.K, res.LogVP_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.VP=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    
    %Predict WS values
    %case {'ws','logws'}
    [Lia,Locb] =ismember({'ws','logws'},lower(prop));
    if find(Lia)
        
        
        
        Desc=train.WS.Desc;
        
        
        if verbose>0
            disp('Predicting LogWS values (Log10 M)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogWS values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.WS.Desc_i);
        
        pred = nnrpred(Xtest,train.WS.model.set.train,train.WS.model.set.y,train.WS.model.set.K,train.WS.model.set.dist_type,train.WS.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogWS_exp=NaN(size(Xtest,1),1);
        end
        res.LogWS_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.WS.model.set.train,Xtest,'auto');
        res.AD_WS=abs(AD.inorout-1)';
        res.AD_WS(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_WS=zeros(size(Xtest,1),1);
        res.Conf_index_WS=zeros(size(Xtest,1),1);
        
        LogWS_CAS_neighbor=cell(size(Xtest,1),5);
        LogWS_InChiKey_neighbor=cell(size(Xtest,1),5);
        LogWS_DTXSID_neighbor=cell(size(Xtest,1),5);
        LogWS_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogWS_Exp_neighbor=zeros(size(Xtest,1),5);
        LogWS_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.WS.CAS);
                    if Li
                        res.LogWS_exp(i,1)=train.WS.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.WS.DTXSID);
                        if Li
                            res.LogWS_exp(i)=train.WS.model.set.y(Lo);
                        end
                    end
                end
            end
            LogWS_CAS_neighbor(i,:)=train.WS.CAS(pred.neighbors(i,:));
            LogWS_InChiKey_neighbor(i,:)=train.WS.InChiKey(pred.neighbors(i,:));
            LogWS_DTXSID_neighbor(i,:)=train.WS.DTXSID(pred.neighbors(i,:));
            LogWS_DSSTOXMPID_neighbor(i,:)=train.WS.DSSTOXMPID(pred.neighbors(i,:));
            LogWS_Exp_neighbor(i,:)=train.WS.model.set.y(pred.neighbors(i,:));
            LogWS_pred_neighbor(i,:)=train.WS.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogWS_Exp_neighbor(i,:),res.LogWS_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_WS(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_WS(i,1)=((1/(1+sqrt(((LogWS_Exp_neighbor(i,:)-LogWS_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_WS(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogWS_pred(i,1)=NaN;
                res.AD_WS(i)=0;
                res.AD_index_WS(i)=0;
                res.Conf_index_WS(i,1)=0;
            end
            if neighbors==1
                res.LogWS_CAS_neighbor(i,:)=LogWS_CAS_neighbor(i,:);
                res.LogWS_InChiKey_neighbor(i,:)=LogWS_InChiKey_neighbor(i,:);
                res.LogWS_DTXSID_neighbor(i,:)=LogWS_DTXSID_neighbor(i,:);
                res.LogWS_DSSTOXMPID_neighbor(i,:)=LogWS_DSSTOXMPID_neighbor(i,:);
                res.LogWS_Exp_neighbor(i,:)=LogWS_Exp_neighbor(i,:);
                res.LogWS_pred_neighbor(i,:)=LogWS_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogWS experimental= %.3f\n', res.LogWS_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogWS predicted= %.3f\n', res.LogWS_pred(i));
                if res.AD_WS(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_WS(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_WS(i));
                %CAS=strjoin(res.LogWS_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.WS.model.set.K, res.LogWS_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.WS.model.set.K, res.LogWS_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.WS.model.set.K, res.LogWS_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogWS experimental= %.3f\n', res.LogWS_exp(i));
                end
                fprintf(output,'LogWS predicted= %.3f\n', res.LogWS_pred(i));
                if res.AD_WS(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_WS(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_WS(i));
                %CAS=strjoin(res.LogWS_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.WS.model.set.K, res.LogWS_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.WS.model.set.K, res.LogWS_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.WS.model.set.K, res.LogWS_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.WS=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    
    %Predict HL values
    %case {'hl','loghl'}
    [Lia,Locb] =ismember({'hl','loghl'},lower(prop));
    if find(Lia)
        
        
        Desc=train.HL.Desc;
        
        
        if verbose>0
            disp('Predicting LogHL values (Log10 atm-m3/V)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogHL values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.HL.Desc_i);
        
        pred = nnrpred(Xtest,train.HL.model.set.train,train.HL.model.set.y,train.HL.model.set.K,train.HL.model.set.dist_type,train.HL.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogHL_exp=NaN(size(Xtest,1),1);
        end
        res.LogHL_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.HL.model.set.train,Xtest,'auto');
        res.AD_HL=abs(AD.inorout-1)';
        res.AD_HL(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_HL=zeros(size(Xtest,1),1);
        res.Conf_index_HL=zeros(size(Xtest,1),1);
        
        HL_CAS_neighbor=cell(size(Xtest,1),5);
        HL_InChiKey_neighbor=cell(size(Xtest,1),5);
        HL_DTXSID_neighbor=cell(size(Xtest,1),5);
        HL_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogHL_Exp_neighbor=zeros(size(Xtest,1),5);
        LogHL_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.HL.CAS);
                    if Li
                        res.LogHL_exp(i,1)=train.LogP.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.HL.DTXSID);
                        if Li
                            res.LogHL_exp(i)=train.HL.model.set.y(Lo);
                        end
                    end
                end
            end
            HL_CAS_neighbor(i,:)=train.HL.CAS(pred.neighbors(i,:));
            HL_InChiKey_neighbor(i,:)=train.HL.InChiKey(pred.neighbors(i,:));
            HL_DTXSID_neighbor(i,:)=train.HL.DTXSID(pred.neighbors(i,:));
            HL_DSSTOXMPID_neighbor(i,:)=train.HL.DSSTOXMPID(pred.neighbors(i,:));
            LogHL_Exp_neighbor(i,:)=train.HL.model.set.y(pred.neighbors(i,:));
            LogHL_pred_neighbor(i,:)=train.HL.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogHL_Exp_neighbor(i,:),res.LogHL_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_HL(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_HL(i,1)=((1/(1+sqrt(((LogHL_Exp_neighbor(i,:)-LogHL_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_HL(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogHL_pred(i,1)=NaN;
                res.AD_HL(i)=0;
                res.AD_index_HL(i)=0;
                res.Conf_index_HL(i,1)=0;
            end
            if neighbors==1
                res.HL_CAS_neighbor(i,:)=HL_CAS_neighbor(i,:);
                res.HL_InChiKey_neighbor(i,:)=HL_InChiKey_neighbor(i,:);
                res.HL_DTXSID_neighbor(i,:)=HL_DTXSID_neighbor(i,:);
                res.HL_DSSTOXMPID_neighbor(i,:)=HL_DSSTOXMPID_neighbor(i,:);
                res.LogHL_Exp_neighbor(i,:)=LogHL_Exp_neighbor(i,:);
                res.LogHL_pred_neighbor(i,:)=LogHL_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogHL experimental= %.3f\n', res.LogHL_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogHL predicted= %.3f\n', res.LogHL_pred(i));
                if res.AD_HL(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_HL(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_HL(i));
                %CAS=strjoin(res.HL_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.HL.model.set.K, res.HL_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.HL.model.set.K, res.LogHL_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.HL.model.set.K, res.LogHL_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogHL experimental= %.3f\n', res.LogHL_exp(i));
                end
                fprintf(output,'LogHL predicted= %.3f\n', res.LogHL_pred(i));
                if res.AD_HL(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_HL(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_HL(i));
                %CAS=strjoin(res.HL_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.HL.model.set.K, res.HL_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.HL.model.set.K, res.LogHL_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.HL.model.set.K, res.LogHL_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.HL=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    
    %Predict RT values
    %case {'rt'}
    [Lia,Locb] =ismember('rt',lower(prop));
    if find(Lia)
        
        
        Desc=train.RT.Desc;
        
        
        if verbose>0
            disp('Predicting RT values (Mins.)...');
            if verbose>1
                disp(['PLS model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting RT values... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.RT.Desc_i);
        
        pred = nnrpred(Xtest,train.RT.model.set.train,train.RT.model.set.y,train.RT.model.set.K,train.RT.model.set.dist_type,train.RT.model.set.scal);
        pred.D=diag(pred.D);
        predpls=plstest(Xtest,train.RT.model);
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.RT_exp=NaN(size(Xtest,1),1);
        end
        res.RT_pred(:,1)=predpls.yc;
        AD=classical_leverage(train.RT.model.set.train,Xtest,'auto');
        res.AD_RT=abs(AD.inorout-1)';
        res.AD_RT(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_RT=zeros(size(Xtest,1),1);
        res.Conf_index_RT=zeros(size(Xtest,1),1);
        
        RT_CAS_neighbor=cell(size(Xtest,1),5);
        RT_DTXSID_neighbor=cell(size(Xtest,1),5);
        RT_Exp_neighbor=zeros(size(Xtest,1),5);
        RT_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.RT.CAS);
                    if Li
                        res.RT_exp(i,1)=train.RT.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.RT.DTXSID);
                        if Li
                            res.RT_exp(i)=train.RT.model.set.y(Lo);
                        end
                    end
                end
            end
            RT_CAS_neighbor(i,:)=train.RT.CAS(pred.neighbors(i,:));
            RT_DTXSID_neighbor(i,:)=train.RT.DTXSID(pred.neighbors(i,:));
            RT_Exp_neighbor(i,:)=train.RT.model.set.y(pred.neighbors(i,:));
            RT_pred_neighbor(i,:)=train.RT.model.yc(pred.neighbors(i,:));
            
            
            res.AD_index_RT(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_RT(i,1)=((1/(1+sqrt(((RT_Exp_neighbor(i,:)-RT_pred_neighbor(i,:)).^2)*pred.w(i,:)')/4.5))+res.AD_index_RT(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.RT_pred(i,1)=NaN;
                res.AD_RT(i)=0;
                res.AD_index_RT(i)=0;
                res.Conf_index_RT(i,1)=0;
            end
            if res.RT_pred(i,1)<0
                res.RT_pred(i,1)=0;
                res.AD_RT(i)=0;
            end
            if neighbors==1
                res.RT_CAS_neighbor(i,:)=RT_CAS_neighbor(i,:);
                res.RT_DTXSID_neighbor(i,:)=RT_DTXSID_neighbor(i,:);
                res.RT_Exp_neighbor(i,:)=RT_Exp_neighbor(i,:);
                res.RT_pred_neighbor(i,:)=RT_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'RT experimental= %.3f\n', res.RT_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'RT predicted= %.3f\n', res.RT_pred(i));
                if res.AD_RT(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_RT(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_RT(i));
                %CAS=strjoin(res.RT_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.RT.model.set.K, res.RT_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.RT.model.set.K, res.RT_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.RT.model.set.K, res.RT_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'RT experimental= %.3f\n', res.RT_exp(i));
                end
                fprintf(output,'RT predicted= %.3f\n', res.RT_pred(i));
                if res.AD_RT(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_RT(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_RT(i));
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.RT.model.set.K, res.RT_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.RT.model.set.K, res.RT_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.RT.model.set.K, res.RT_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.RT=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('predpls');
        clear('AD');
        %end clean memory
        
    end
    
    %Predict KOA values
    %case {'koa','logkoa'}
    [Lia,Locb] =ismember({'koa','logkoa'},lower(prop));
    if find(Lia)
        
        Desc=train.KOA.Desc;
        
        
        if verbose>0
            disp('Predicting LogKOA values (Log10)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogKOA values... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.KOA.Desc_i);
        
        pred = nnrpred(Xtest,train.KOA.model.set.train,train.KOA.model.set.y,train.KOA.model.set.K,train.KOA.model.set.dist_type,train.KOA.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogKOA_exp=NaN(size(Xtest,1),1);
        end
        res.LogKOA_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.KOA.model.set.train,Xtest,'auto');
        res.AD_KOA=abs(AD.inorout-1)';
        res.AD_KOA(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_KOA=zeros(size(Xtest,1),1);
        res.Conf_index_KOA=zeros(size(Xtest,1),1);
        
        KOA_CAS_neighbor=cell(size(Xtest,1),5);
        KOA_InChiKey_neighbor=cell(size(Xtest,1),5);
        KOA_DTXSID_neighbor=cell(size(Xtest,1),5);
        KOA_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogKOA_Exp_neighbor=zeros(size(Xtest,1),5);
        LogKOA_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.KOA.CAS);
                    if Li
                        res.LogKOA_exp(i,1)=train.KOA.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.KOA.DTXSID);
                        if Li
                            res.LogKOA_exp(i)=train.KOA.model.set.y(Lo);
                        end
                    end
                end
            end
            KOA_CAS_neighbor(i,:)=train.KOA.CAS(pred.neighbors(i,:));
            KOA_InChiKey_neighbor(i,:)=train.KOA.InChiKey(pred.neighbors(i,:));
            KOA_DTXSID_neighbor(i,:)=train.KOA.DTXSID(pred.neighbors(i,:));
            KOA_DSSTOXMPID_neighbor(i,:)=train.KOA.DSSTOXMPID(pred.neighbors(i,:));
            LogKOA_Exp_neighbor(i,:)=train.KOA.model.set.y(pred.neighbors(i,:));
            LogKOA_pred_neighbor(i,:)=train.KOA.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogKOA_Exp_neighbor(i,:),res.LogKOA_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_KOA(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_KOA(i,1)=((1/(1+sqrt(((LogKOA_Exp_neighbor(i,:)-LogKOA_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_KOA(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogKOA_pred(i,1)=NaN;
                res.AD_KOA(i)=0;
                res.AD_index_KOA(i)=0;
                res.Conf_index_KOA(i,1)=0;
            end
            if neighbors==1
                res.KOA_CAS_neighbor(i,:)=KOA_CAS_neighbor(i,:);
                res.KOA_InChiKey_neighbor(i,:)=KOA_InChiKey_neighbor(i,:);
                res.KOA_DTXSID_neighbor(i,:)=KOA_DTXSID_neighbor(i,:);
                res.KOA_DSSTOXMPID_neighbor(i,:)=KOA_DSSTOXMPID_neighbor(i,:);
                res.LogKOA_Exp_neighbor(i,:)=LogKOA_Exp_neighbor(i,:);
                res.LogKOA_pred_neighbor(i,:)=LogKOA_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogKOA experimental= %.3f\n', res.LogKOA_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogKOA predicted= %.3f\n', res.LogKOA_pred(i));
                if res.AD_KOA(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_KOA(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_KOA(i));
                %CAS=strjoin(res.KOA_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KOA.model.set.K, res.KOA_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KOA.model.set.K, res.LogKOA_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KOA.model.set.K, res.LogKOA_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogKOA experimental= %.3f\n', res.LogKOA_exp(i));
                end
                fprintf(output,'LogKOA predicted= %.3f\n', res.LogKOA_pred(i));
                if res.AD_KOA(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_KOA(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_KOA(i));
                %CAS=strjoin(res.KOA_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KOA.model.set.K, res.KOA_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KOA.model.set.K, res.LogKOA_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KOA.model.set.K, res.LogKOA_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.KOA=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    %Predict pka values
    %case {'pka'}
    [Lia,Locb] =ismember({'pka','logd'},lower(prop));
    if find(Lia)
        
        Desc=train.pKa.Desc;
        %             Desc_a=train.pKa.Desc_a;
        %             Desc_b=train.pKa.Desc_b;
        
        if verbose>0
            disp('Predicting pKa values (unitless)...');
            if verbose>1
                disp(['SVM models with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0 && Lia(1)
            fprintf(output,'\n\n\t\t\t\t\t Predicting pKa values... \n\n			==============================================================  \n\n');
        end
        
        
        Xtest=Xin(:,train.pKa.Desc_i);
        Xtest_a=table2array(XinFP(:,train.pKa.Desc_ai));
        Xtest_b=table2array(XinFP(:,train.pKa.Desc_bi));
        
        
        pred = knnpred(Xtest,train.pKa.model.set.train,train.pKa.model.set.class,train.pKa.model.set.K,train.pKa.model.set.dist_type,train.pKa.model.set.param.pret_type);
        pred.D=diag(pred.D);
        pKa_a(:,1)=svmpredict([1:1:length(Xtest_a(:,1))]',Xtest_a,train.pKa.model_a,'-q');
        %AD_a = nnrpred(Xtest_a,train.pKa_a.model.set.train,train.pKa_a.model.set.y,train.pka_a.model.set.K,train.pka_a.model.set.dist_type,train.pka_a.model.set.param.pret_type);
        
        pKa_b(:,1)=svmpredict([1:1:length(Xtest_b(:,1))]',Xtest_b,train.pKa.model_b,'-q');
        %AD_b = nnrpred(Xtest_b,train.pka_b.model.set.train,train.pka_b.model.set.y,train.pka_b.model.set.K,train.pka_b.model.set.dist_type,train.pka_b.model.set.param.pret_type);
        
        res.MoleculeID=MoleculeNames;
        res.ionization=zeros(size(Xtest,1),1);
        pKa_ac_ba_amp=pred.class_pred;
        res.pKa_a_pred=pKa_a;
        res.pKa_b_pred=pKa_b;
        
        AD=classical_leverage(train.pKa.model.set.train,Xtest,'auto');
        res.AD_pKa=abs(AD.inorout-1)';
        res.AD_pKa(round(pred.D,3)==0)=1;
        
        
        res.AD_index_pKa=zeros(size(Xtest,1),1);
        res.Conf_index_pKa=zeros(size(Xtest,1),1);
        
        pKa_CAS_neighbor=cell(size(Xtest,1),3);
        pKa_InChiKey_neighbor=cell(size(Xtest,1),3);
        pKa_DTXSID_neighbor=cell(size(Xtest,1),3);
        %pKa_DSSTOXMPID_neighbor=cell(size(Xtest,1),3);
        pKa_Exp_neighbor=zeros(size(Xtest,1),3);
        pKa_pred_neighbor=zeros(size(Xtest,1),3);
        
        
        for i=1:size(Xtest(:,1))
            if sum(Xin(i,13:14))-sum(Xin(i,731:732))==0
                pKa_ac_ba_amp(i)=NaN;
                res.ionization(i)=0;
                res.pKa_a_pred(i)=NaN;
                res.pKa_b_pred(i)=NaN;
                
            else
                
                if pred.class_pred(i)==1
                    res.pKa_b_pred(i,1)=NaN;
                    res.ionization(i)=1;
                elseif pred.class_pred(i)==2
                    res.pKa_a_pred(i,1)=NaN;
                    res.ionization(i)=1;
                elseif pred.class_pred(i)==3
                    res.ionization(i)=2;
                    
                end
            end
            
            
            pKa_CAS_neighbor(i,:)=train.pKa.CAS(pred.neighbors(i,:));
            pKa_InChiKey_neighbor(i,:)=train.pKa.InChiKey(pred.neighbors(i,:));
            pKa_DTXSID_neighbor(i,:)=train.pKa.DTXSID(pred.neighbors(i,:));
            %                 pKa_DSSTOXMPID_neighbor(i,:)=train.pKa.DSSTOXMPID(pred.neighbors(i,:));
            pKa_Exp_neighbor(i,:)=train.pKa.model.set.y(pred.neighbors(i,:));
            pKa_pred_neighbor(i,:)=train.pKa.model.set.yc(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogKOA_Exp_neighbor(i,:),res.LogKOA_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_pKa(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            res.Conf_index_pKa(i,1)=((1/(1+sqrt(((pKa_Exp_neighbor(i,:)-pKa_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_pKa(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.pKa_a_pred(i,1)=NaN;
                res.pKa_b_pred(i,1)=NaN;
                res.AD_pKa(i,1)=0;
                res.AD_index_pKa(i,1)=0;
                res.Conf_index_pKa(i,1)=0;
            end
            if neighbors==1
                res.pKa_CAS_neighbor(i,:)=pKa_CAS_neighbor(i,:);
                res.pKa_InChiKey_neighbor(i,:)=pKa_InChiKey_neighbor(i,:);
                res.pKa_DTXSID_neighbor(i,:)=pKa_DTXSID_neighbor(i,:);
                %res.pKa_DSSTOXMPID_neighbor(i,:)=pKa_DSSTOXMPID_neighbor(i,:);
                res.pKa_Exp_neighbor(i,:)=pKa_Exp_neighbor(i,:);
                res.pKa_pred_neighbor(i,:)=pKa_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1 && Lia(1)
                %res.Xtest=Xtest;
                fprintf(output(Locb(1)),'\t Molecule %s:\n', MoleculeNames{i});
%                 if exp
%                     fprintf(output(Locb(1)),'BP experimental= %.3f\n', res.BP_exp(i));
%                 end
                fprintf(output(Locb(1)),'pKa acidic and basic predicted= %.3f, %.3f\n', res.pKa_a_pred(i),res.pKa_b_pred(i));
                if res.AD_pKa(i)==1
                    fprintf(output(Locb(1)),'AD: inside\n');
                else
                    fprintf(output(Locb(1)),'AD: outside\n');
                end
                fprintf(output(Locb(1)),'AD_index= %.2f\n', res.AD_index_pKa(i));
                fprintf(output(Locb(1)),'Conf_index= %.2f\n', res.Conf_index_pKa(i));
                %CAS=strjoin(res.RT_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(1)),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.pKa.model.set.K, res.pKa_CAS_neighbor{i,1:3});
                    fprintf(output(Locb(1)),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.pKa.model.set.K, res.pKa_Exp_neighbor(i,1:3));
                    fprintf(output(Locb(1)),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.pKa.model.set.K, res.pKa_pred_neighbor(i,1:3));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0 && Lia(1)
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                fprintf(output,'pKa acidic and basic predicted= %.3f, %.3f\n', res.pKa_a_pred(i),res.pKa_b_pred(i));
                if res.AD_pKa(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_pKa(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_pKa(i));
                %CAS=strjoin(res.KOA_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.pKa.model.set.K, res.pKa_CAS_neighbor{i,1:3});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.pKa.model.set.K, res.pKa_Exp_neighbor(i,1:3));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.pKa.model.set.K, res.pKa_pred_neighbor(i,1:3));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv') && Lia(1)
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(1)},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(1)));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv') && Lia(1)
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.pKa=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('XinFP');
        clear('Xtest_a');
        clear('Xtest_b');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    %Predict LogD values
    %case 'logd'
    [Lia,Locb] =ismember('logd',lower(prop));
    if find(Lia)
        
        
        if verbose>0
            disp('Predicting LogD values (Log10)...');
            if verbose>1
                disp('Predictions based on pKa and LogP ');
            end
        end
        res.MoleculeID=MoleculeNames;
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogD values... \n\n			==============================================================  \n\n');
        end
        if sep==1
            
            res.LogD55_pred=resf.LogP.LogP_pred;
            res.LogD74_pred=resf.LogP.LogP_pred;
            res.AD_LogD=resf.LogP.AD_LogP+resf.pKa.AD_pKa;
            res.AD_LogD(find(res.AD_LogD==1))=0;
            res.AD_LogD(find(res.AD_LogD==2))=1;
            res.AD_index_LogD=0.5*resf.pKa.AD_index_pKa+0.5*resf.LogP.AD_index_LogP;
            res.Conf_index_LogD=0.5*resf.pKa.Conf_index_pKa+0.5*resf.LogP.Conf_index_LogP;
            
            if neighbors==1
                res.LogD_CAS_neighbor=resf.LogP.LogP_CAS_neighbor;
                res.LogD_InChiKey_neighbor=resf.LogP.LogP_InChiKey_neighbor;
                res.LogD_DTXSID_neighbor=resf.LogP.LogP_DTXSID_neighbor;
                res.LogD_DSSTOXMPID_neighbor=resf.LogP.LogP_DSSTOXMPID_neighbor;
                %res.LogD_Exp_neighbor=LogP_Exp_neighbor;
                %res.LogD_pred_neighbor=LogP_pred_neighbor;
            end
            
            for i=1:length(res.LogD55_pred)
                
                if pKa_ac_ba_amp(i)==1
                    res.LogD55_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^(5.5-resf.pKa.pKa_a_pred(i,1)));
                    res.LogD74_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^(7.4-resf.pKa.pKa_a_pred(i,1)));
                elseif pKa_ac_ba_amp(i)==2
                    res.LogD55_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^(resf.pKa.pKa_b_pred(i,1)-5.5));
                    res.LogD74_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^(resf.pKa.pKa_b_pred(i,1)-7.4));
                elseif pKa_ac_ba_amp(i)==3
                    res.LogD55_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^abs(0.5*resf.pKa.pKa_a_pred(i,1)+0.5*resf.pKa.pKa_b_pred(i,1)-5.5));
                    res.LogD74_pred(i,1)=resf.LogP.LogP_pred(i,1)-log10(1+10^abs(0.5*resf.pKa.pKa_a_pred(i,1)+0.5*resf.pKa.pKa_b_pred(i,1)-7.4));
                end
                
                if strcmpi(ext,'.txt')
                    fprintf(output(Locb(1)),'\t Molecule %s:\n', res.MoleculeID{i});
                    fprintf(output(Locb(1)),'LogD pH 5.5 predicted= %.3f\n', res.LogD55_pred(i));
                    fprintf(output(Locb(1)),'LogD pH 7.4 predicted= %.3f\n', res.LogD74_pred(i));
                    if res.AD_LogD(i)==1
                        fprintf(output(Locb(1)),'AD: inside\n');
                    else
                        fprintf(output(Locb(1)),'AD: outside\n');
                    end
                    fprintf(output(Locb(1)),'AD_index= %.2f\n', res.AD_index_LogD(i));
                    fprintf(output(Locb(1)),'Conf_index= %.2f\n', res.Conf_index_LogD(i));
                    if neighbors==1
                        fprintf(output(Locb(1)),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.LogP.model.set.K,resf.LogP.LogP_CAS_neighbor{i,1:5});
                        %fprintf(output(Locb(1)),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.LogP.model.set.K, res.LogP_Exp_neighbor(i,1:5));
                        %fprintf(output(Locb(1)),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.LogP.model.set.K, res.LogP_pred_neighbor(i,1:5));
                    end
                end
            end
            
            if strcmpi(ext,'.csv')
                T=struct2table(res);
                writetable(T,FileOut{Locb(1)},'Delimiter',',');%,'QuoteStrings',true);
                fclose(output(Locb(1)));
                clear('T');
            end
            resf.LogD=res;
            clear('res');
            
            
        else
            res.LogD55_pred=res.LogP_pred;
            res.LogD74_pred=res.LogP_pred;
            res.AD_LogD=res.AD_LogP+res.AD_pKa;
            res.AD_LogD(find(res.AD_LogD==1))=0;
            res.AD_LogD(find(res.AD_LogD==2))=1;
            res.AD_index_LogD=0.5*res.AD_index_pKa+0.5*res.AD_index_LogP;
            res.Conf_index_LogD=0.5*res.Conf_index_pKa+0.5*res.Conf_index_LogP;
            
            if neighbors==1
                res.LogD_CAS_neighbor=res.LogP_CAS_neighbor;
                res.LogD_InChiKey_neighbor=res.LogP_InChiKey_neighbor;
                res.LogD_DTXSID_neighbor=res.LogP_DTXSID_neighbor;
                res.LogD_DSSTOXMPID_neighbor=res.LogP_DSSTOXMPID_neighbor;
                %res.LogD_Exp_neighbor=LogP_Exp_neighbor;
                %res.LogD_pred_neighbor=LogP_pred_neighbor;
            end
            
            for i=1:length(res.LogD55_pred)
                if pKa_ac_ba_amp(i)==1
                    res.LogD55_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^(5.5-res.pKa_a_pred(i,1)));
                    res.LogD74_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^(7.4-res.pKa_a_pred(i,1)));
                elseif pKa_ac_ba_amp(i)==2
                    res.LogD55_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^(res.pKa_b_pred(i,1)-5.5));
                    res.LogD74_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^(res.pKa_b_pred(i,1)-7.4));
                elseif pKa_ac_ba_amp(i)==3
                    res.LogD55_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^abs(0.5*res.pKa_a_pred(i,1)+0.5*res.pKa_b_pred(i,1)-5.5));
                    res.LogD74_pred(i,1)=res.LogP_pred(i,1)-log10(1+10^abs(0.5*res.pKa_a_pred(i,1)+0.5*res.pKa_b_pred(i,1)-7.4));
                end
                
                if strcmpi(ext,'.txt')
                    fprintf(output,'\t Molecule %s:\n',res.MoleculeID{i});
                    fprintf(output,'LogD pH 5.5 predicted= %.3f\n', res.LogD55_pred(i));
                    fprintf(output,'LogD pH 7.4 predicted= %.3f\n', res.LogD74_pred(i));
                    if res.AD_LogD(i)==1
                        fprintf(output,'AD: inside\n');
                    else
                        fprintf(output,'AD: outside\n');
                    end
                    fprintf(output,'AD_index= %.2f\n', res.AD_index_LogD(i));
                    fprintf(output,'Conf_index= %.2f\n', res.Conf_index_LogD(i));
                    if neighbors==1
                        fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.LogP.model.set.K, res.LogP_CAS_neighbor{i,1:5});
                        %fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.LogP.model.set.K, res.LogP_Exp_neighbor(i,1:5));
                        %fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.LogP.model.set.K, res.LogP_pred_neighbor(i,1:5));
                    end
                    
                end
                
                
            end
        end
        
    end
    
    %  Env. Fate Endpoints
    
    if verbose> 0 && (ef||all)
        fprintf(1,'---------- Env. Fate Endpoints ----------\n');
    end
    %Predict AOH values
    %case {'aop','logoh','aoh'}
    [Lia,Locb] =ismember({'aop','logoh','aoh'},lower(prop));
    if find(Lia)
        
        Desc=train.AOH.Desc;
        
        
        if verbose>0
            disp('Predicting LogOH values (Log10 cm3/molecule-sec)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogOH values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.AOH.Desc_i);
        
        pred = nnrpred(Xtest,train.AOH.model.set.train,train.AOH.model.set.y,train.AOH.model.set.K,train.AOH.model.set.dist_type,train.AOH.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogOH_exp=NaN(size(Xtest,1),1);
        end
        res.LogOH_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.AOH.model.set.train,Xtest,'auto');
        res.AD_AOH=abs(AD.inorout-1)';
        res.AD_AOH(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_AOH=zeros(size(Xtest,1),1);
        res.Conf_index_AOH=zeros(size(Xtest,1),1);
        
        AOH_CAS_neighbor=cell(size(Xtest,1),5);
        AOH_InChiKey_neighbor=cell(size(Xtest,1),5);
        AOH_DTXSID_neighbor=cell(size(Xtest,1),5);
        AOH_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogOH_Exp_neighbor=zeros(size(Xtest,1),5);
        LogOH_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.AOH.CAS);
                    if Li
                        res.LogOH_exp(i,1)=train.AOH.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.AOH.DTXSID);
                        if Li
                            res.LogOH_exp(i)=train.AOH.model.set.y(Lo);
                        end
                    end
                end
            end
            AOH_CAS_neighbor(i,:)=train.AOH.CAS(pred.neighbors(i,:));
            AOH_InChiKey_neighbor(i,:)=train.AOH.InChiKey(pred.neighbors(i,:));
            AOH_DTXSID_neighbor(i,:)=train.AOH.DTXSID(pred.neighbors(i,:));
            AOH_DSSTOXMPID_neighbor(i,:)=train.AOH.DSSTOXMPID(pred.neighbors(i,:));
            LogOH_Exp_neighbor(i,:)=train.AOH.model.set.y(pred.neighbors(i,:));
            LogOH_pred_neighbor(i,:)=train.AOH.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogOH_Exp_neighbor(i,:),res.LogOH_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_AOH(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_AOH(i,1)=((1/(1+sqrt(((LogOH_Exp_neighbor(i,:)-LogOH_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_AOH(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogOH_pred(i,1)=NaN;
                res.AD_AOH(i)=0;
                res.AD_index_AOH(i)=0;
                res.Conf_index_AOH(i,1)=0;
            end
            if neighbors==1
                res.AOH_CAS_neighbor(i,:)=AOH_CAS_neighbor(i,:);
                res.AOH_InChiKey_neighbor(i,:)=AOH_InChiKey_neighbor(i,:);
                res.AOH_DTXSID_neighbor(i,:)=AOH_DTXSID_neighbor(i,:);
                res.AOH_DSSTOXMPID_neighbor(i,:)=AOH_DSSTOXMPID_neighbor(i,:);
                res.LogOH_Exp_neighbor(i,:)=LogOH_Exp_neighbor(i,:);
                res.LogOH_pred_neighbor(i,:)=LogOH_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogOH experimental= %.3f\n', res.LogOH_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogOH predicted= %.3f\n', res.LogOH_pred(i));
                if res.AD_AOH(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_AOH(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_AOH(i));
                %CAS=strjoin(res.AOH_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.AOH.model.set.K, res.AOH_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.AOH.model.set.K, res.LogOH_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.AOH.model.set.K, res.LogOH_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0

                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogOH experimental= %.3f\n', res.LogOH_exp(i));
                end
                fprintf(output,'LogOH predicted= %.3f\n', res.LogOH_pred(i));
                if res.AD_AOH(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_AOH(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_AOH(i));
                %CAS=strjoin(res.AOH_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.AOH.model.set.K, res.AOH_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.AOH.model.set.K, res.LogOH_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.AOH.model.set.K, res.LogOH_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.AOH=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    
    %Predict BCF values
    %case {'bcf', 'logbcf'}
    [Lia,Locb] =ismember({'bcf','logbcf'},lower(prop));
    if find(Lia)
        
        
        
        Desc=train.BCF.Desc;
        
        
        if verbose>0
            disp('Predicting LogBCF values (Log10)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogBCF values... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.BCF.Desc_i);
        
        pred = nnrpred(Xtest,train.BCF.model.set.train,train.BCF.model.set.y,train.BCF.model.set.K,train.BCF.model.set.dist_type,train.BCF.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogBCF_exp=NaN(size(Xtest,1),1);
        end
        res.LogBCF_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.BCF.model.set.train,Xtest,'auto');
        res.AD_BCF=abs(AD.inorout-1)';
        res.AD_BCF(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_BCF=zeros(size(Xtest,1),1);
        res.Conf_index_BCF=zeros(size(Xtest,1),1);
        
        LogBCF_CAS_neighbor=cell(size(Xtest,1),5);
        LogBCF_InChiKey_neighbor=cell(size(Xtest,1),5);
        LogBCF_DTXSID_neighbor=cell(size(Xtest,1),5);
        LogBCF_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogBCF_Exp_neighbor=zeros(size(Xtest,1),5);
        LogBCF_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.BCF.CAS);
                    if Li
                        res.LogBCF_exp(i,1)=train.BCF.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.BCF.DTXSID);
                        if Li
                            res.LogBCF_exp(i)=train.BCF.model.set.y(Lo);
                        end
                    end
                end
            end
            LogBCF_CAS_neighbor(i,:)=train.BCF.CAS(pred.neighbors(i,:));
            LogBCF_InChiKey_neighbor(i,:)=train.BCF.InChiKey(pred.neighbors(i,:));
            LogBCF_DTXSID_neighbor(i,:)=train.BCF.DTXSID(pred.neighbors(i,:));
            LogBCF_DSSTOXMPID_neighbor(i,:)=train.BCF.DSSTOXMPID(pred.neighbors(i,:));
            LogBCF_Exp_neighbor(i,:)=train.BCF.model.set.y(pred.neighbors(i,:));
            LogBCF_pred_neighbor(i,:)=train.BCF.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogBCF_Exp_neighbor(i,:),res.LogBCF_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            %res.Conf_index2(i,1)=(res.Conf_index(i)*res.AD_index(i))^0.5;
            
            res.AD_index_BCF(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_BCF(i,1)=((1/(1+sqrt(((LogBCF_Exp_neighbor(i,:)-LogBCF_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_BCF(i,1))/2;
            if isnan(res.AD_index_BCF(i))
                res.LogBCF_pred(i,1)=NaN;
                res.AD_BCF(i)=0;
                res.AD_index_BCF(i)=0;
                res.Conf_index_BCF(i,1)=0;
            end
            
            
            if neighbors==1
                res.LogBCF_CAS_neighbor(i,:)=LogBCF_CAS_neighbor(i,:);
                res.LogBCF_InChiKey_neighbor(i,:)=LogBCF_InChiKey_neighbor(i,:);
                res.LogBCF_DTXSID_neighbor(i,:)=LogBCF_DTXSID_neighbor(i,:);
                res.LogBCF_DSSTOXMPID_neighbor(i,:)=LogBCF_DSSTOXMPID_neighbor(i,:);
                res.LogBCF_Exp_neighbor(i,:)=LogBCF_Exp_neighbor(i,:);
                res.LogBCF_pred_neighbor(i,:)=LogBCF_pred_neighbor(i,:);
            end
            
            
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogBCF experimental= %.3f\n', res.LogBCF_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogBCF predicted= %.3f\n', res.LogBCF_pred(i));
                if res.AD_BCF(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_BCF(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_BCF(i));
                %CAS=strjoin(res.LogBCF_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.BCF.model.set.K, res.LogBCF_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.BCF.model.set.K, res.LogBCF_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.BCF.model.set.K, res.LogBCF_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogBCF experimental= %.3f\n', res.LogBCF_exp(i));
                end
                fprintf(output,'LogBCF predicted= %.3f\n', res.LogBCF_pred(i));
                if res.AD_BCF(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_BCF(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_BCF(i));
                %CAS=strjoin(res.LogBCF_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.BCF.model.set.K, res.LogBCF_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.BCF.model.set.K, res.LogBCF_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.BCF.model.set.K, res.LogBCF_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.BCF=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    %Predict Biodegradability values
    %case {'biohc','biohl','biodeg','biodeghl'}
    [Lia,Locb] =ismember({'biohc','biohl','biodeg','biodeghl'},lower(prop));
    if find(Lia)
        
        Desc=train.Biodeg.Desc;
        
        
        if verbose>0
            disp('Predicting Biodeg. half-life values (Log10 days)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
        end
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting Biodegradability in LogHalfLife... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.Biodeg.Desc_i);
        
        pred = nnrpred(Xtest,train.Biodeg.model.set.train,train.Biodeg.model.set.y,train.Biodeg.model.set.K,train.Biodeg.model.set.dist_type,train.Biodeg.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.BioDeg_exp=NaN(size(Xtest,1),1);
        end
        res.BioDeg_LogHalfLife_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.Biodeg.model.set.train,Xtest,'auto');
        res.AD_BioDeg=abs(AD.inorout-1)';
        res.AD_BioDeg(round(pred.D,3)==0)=1;
        
        
        %             res.dc=pred.dc;
        %res.AD_index1=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index1)
        %                 res.AD_index1=0;
        %             end
        
        res.AD_index_BioDeg=zeros(size(Xtest,1),1);
        res.Conf_index_BioDeg=zeros(size(Xtest,1),1);
        
        BioDeg_CAS_neighbor=cell(size(Xtest,1),5);
        BioDeg_InChiKey_neighbor=cell(size(Xtest,1),5);
        BioDeg_DTXSID_neighbor=cell(size(Xtest,1),5);
        BioDeg_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        BioDeg_LogHalfLife_Exp_neighbor=zeros(size(Xtest,1),5);
        BioDeg_LogHalfLife_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.Biodeg.CAS);
                    if Li
                        res.BioDeg_exp(i,1)=train.Biodeg.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.Biodeg.DTXSID);
                        if Li
                            res.BioDeg_exp(i)=train.Biodeg.model.set.y(Lo);
                        end
                    end
                end
            end
            BioDeg_CAS_neighbor(i,:)=train.Biodeg.CAS(pred.neighbors(i,:));
            BioDeg_InChiKey_neighbor(i,:)=train.Biodeg.InChiKey(pred.neighbors(i,:));
            BioDeg_DTXSID_neighbor(i,:)=train.Biodeg.DTXSID(pred.neighbors(i,:));
            BioDeg_DSSTOXMPID_neighbor(i,:)=train.Biodeg.DSSTOXMPID(pred.neighbors(i,:));
            BioDeg_LogHalfLife_Exp_neighbor(i,:)=train.Biodeg.model.set.y(pred.neighbors(i,:));
            BioDeg_LogHalfLife_pred_neighbor(i,:)=train.Biodeg.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.BioDeg_LogHalfLife_Exp_neighbor(i,:),res.BioDeg_LogHalfLife_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_BioDeg(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_BioDeg(i,1)=((1/(1+sqrt(((BioDeg_LogHalfLife_Exp_neighbor(i,:)-BioDeg_LogHalfLife_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_BioDeg(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.BioDeg_LogHalfLife_pred(i,1)=NaN;
                res.AD_BioDeg(i)=0;
                res.AD_index_BioDeg(i)=0;
                res.Conf_index_BioDeg(i,1)=0;
            end
            if neighbors==1
                res.BioDeg_CAS_neighbor(i,:)=BioDeg_CAS_neighbor(i,:);
                res.BioDeg_InChiKey_neighbor(i,:)=BioDeg_InChiKey_neighbor(i,:);
                res.BioDeg_DTXSID_neighbor(i,:)=BioDeg_DTXSID_neighbor(i,:);
                res.BioDeg_DSSTOXMPID_neighbor(i,:)=BioDeg_DSSTOXMPID_neighbor(i,:);
                res.BioDeg_LogHalfLife_Exp_neighbor(i,:)=BioDeg_LogHalfLife_Exp_neighbor(i,:);
                res.BioDeg_LogHalfLife_pred_neighbor(i,:)=BioDeg_LogHalfLife_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'BioDeg_LogHalfLife experimental= %.3f\n', res.BioDeg_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'BioDeg_LogHalfLife predicted= %.3f\n', res.BioDeg_LogHalfLife_pred(i));
                if res.AD_BioDeg(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_BioDeg(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_BioDeg(i));
                %CAS=strjoin(res.BioDeg_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.Biodeg.model.set.K, res.BioDeg_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.Biodeg.model.set.K, res.BioDeg_LogHalfLife_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.Biodeg.model.set.K, res.BioDeg_LogHalfLife_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0

                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'BioDeg_LogHalfLife experimental= %.3f\n', res.BioDeg_exp(i));
                end
                fprintf(output,'BioDeg_LogHalfLife predicted= %.3f\n', res.BioDeg_LogHalfLife_pred(i));
                if res.AD_BioDeg(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_BioDeg(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_BioDeg(i));
                %CAS=strjoin(res.BioDeg_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.Biodeg.model.set.K, res.BioDeg_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.Biodeg.model.set.K, res.BioDeg_LogHalfLife_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.Biodeg.model.set.K, res.BioDeg_LogHalfLife_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.BioDeg=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    %Predict RBiodeg values
    %case {'biowin','rb','readybiodeg','rbiodeg'}
    [Lia,Locb] =ismember({'biowin','rb','readybiodeg','rbiodeg'},lower(prop));
    if find(Lia)
        
        Desc=train.RBioDeg.Desc;
        
        if verbose>0
            disp('Predicting Ready-Biodegradability (Binary 0/1)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting Ready Biodegradability... \n\n			============================================================== \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.RBioDeg.Desc_i);
        
        pred = knnpred(Xtest,train.RBioDeg.model.set.train,train.RBioDeg.model.set.class,train.RBioDeg.model.set.K,train.RBioDeg.model.set.dist_type,train.RBioDeg.model.set.param.pret_type);
        pred.D=diag(pred.D);
        %pred.w = (ones(1,train.RBioDeg.model.set.K)./train.RBioDeg.model.set.K)';
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.ReadyBiodeg_exp=NaN(size(Xtest,1),1);
        end
        res.ReadyBiodeg_pred(:,1)=pred.class_pred-1;
        AD=classical_leverage(train.RBioDeg.model.set.train,Xtest,'auto');
        res.AD_ReadyBiodeg=abs(AD.inorout-1)';
        res.AD_ReadyBiodeg(round(pred.D,3)==0)=1;
        %
        
        
        %res.dc=pred.dc;
        res.AD_index_ReadyBiodeg=1./(1+nanmedian(pred.dc,2));
        res.AD_index_ReadyBiodeg(isnan(res.AD_index_ReadyBiodeg))=0;
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        
        
        %             res.AD_index=zeros(size(Xtest,1),1);
        %             res.Conf_index1=zeros(size(Xtest,1),1);
        res.Conf_index_ReadyBiodeg=zeros(size(Xtest,1),1);
        
        ReadyBiodeg_CAS_neighbor=cell(size(Xtest,1),5);
        ReadyBiodeg_InChiKey_neighbor=cell(size(Xtest,1),5);
        ReadyBiodeg_DTXSID_neighbor=cell(size(Xtest,1),5);
        ReadyBiodeg_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        ReadyBiodeg_Exp_neighbor=zeros(size(Xtest,1),5);
        ReadyBiodeg_pred_neighbor=zeros(size(Xtest,1),5);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.RBioDeg.CAS);
                    if Li
                        res.ReadyBiodeg_exp(i,1)=train.RBioDeg.model.set.class(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.RBioDeg.DTXSID);
                        if Li
                            res.ReadyBiodeg_exp(i)=train.RBioDeg.model.set.class(Lo);
                        end
                    end
                end
            end
            ReadyBiodeg_CAS_neighbor(i,:)=train.RBioDeg.CAS(pred.neighbors(i,:));
            ReadyBiodeg_InChiKey_neighbor(i,:)=train.RBioDeg.InChiKey(pred.neighbors(i,:));
            ReadyBiodeg_DTXSID_neighbor(i,:)=train.RBioDeg.DTXSID(pred.neighbors(i,:));
            ReadyBiodeg_DSSTOXMPID_neighbor(i,:)=train.RBioDeg.DSSTOXMPID(pred.neighbors(i,:));
            ReadyBiodeg_Exp_neighbor(i,:)=train.RBioDeg.model.set.class(pred.neighbors(i,:))-1;
            ReadyBiodeg_pred_neighbor(i,:)=train.RBioDeg.model.class_calc(pred.neighbors(i,:))-1;
            
            rmse=calc_reg_param(ReadyBiodeg_Exp_neighbor(i,:),ReadyBiodeg_pred_neighbor(i,:));
            
            
            res.Conf_index_ReadyBiodeg(i,1)=((1/(1+rmse.RMSEC))+res.AD_index_ReadyBiodeg(i))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.ReadyBiodeg_pred(i,1)=NaN;
                res.AD_ReadyBiodeg(i)=0;
                res.AD_index_ReadyBiodeg(i)=0;
                res.Conf_index_ReadyBiodeg(i,1)=0;
            end
            
            %                 res.AD_index(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            %
            %                 if isnan(res.AD_index(i))
            %                     res.AD_index(i)=0;
            %                 end
            
            
            %                res.Conf_index(i,1)=1/(1+sqrt(((res.ReadyBiodeg_Exp_neighbor(i,:)-res.ReadyBiodeg_pred_neighbor(i,:)).^2)*pred.w(i,:)'));
            if neighbors==1
                res.ReadyBiodeg_CAS_neighbor(i,:)=ReadyBiodeg_CAS_neighbor(i,:);
                res.ReadyBiodeg_InChiKey_neighbor(i,:)=ReadyBiodeg_InChiKey_neighbor(i,:);
                res.ReadyBiodeg_DTXSID_neighbor(i,:)=ReadyBiodeg_DTXSID_neighbor(i,:);
                res.ReadyBiodeg_DSSTOXMPID_neighbor(i,:)=ReadyBiodeg_DSSTOXMPID_neighbor(i,:);
                res.ReadyBiodeg_Exp_neighbor(i,:)=ReadyBiodeg_Exp_neighbor(i,:);
                res.ReadyBiodeg_pred_neighbor(i,:)=ReadyBiodeg_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'ReadyBiodeg experimental= %.3f\n', res.ReadyBiodeg_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'ReadyBiodeg predicted= %d\n', res.ReadyBiodeg_pred(i));
                if res.AD_ReadyBiodeg(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_ReadyBiodeg(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_ReadyBiodeg(i));
                %CAS=strjoin(res.ReadyBiodeg_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15d,%15d,%15d,%15d,%15d\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14d,%15d,%15d,%15d,%15d\n\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0

                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'ReadyBiodeg experimental= %.3f\n', res.ReadyBiodeg_exp(i));
                end
                fprintf(output,'ReadyBiodeg predicted= %d\n', res.ReadyBiodeg_pred(i));
                if res.AD_ReadyBiodeg(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_ReadyBiodeg(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_ReadyBiodeg(i));
                %CAS=strjoin(res.ReadyBiodeg_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15d,%15d,%15d,%15d,%15d\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14d,%15d,%15d,%15d,%15d\n\n',train.RBioDeg.model.set.K, res.ReadyBiodeg_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.RBioDeg=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
        
    end
    %Predict KM values
    %case {'km','logkm'}
    [Lia,Locb] =ismember({'km','logkm'},lower(prop));
    if find(Lia)
        
        Desc=train.KM.Desc;
        
        if verbose>0
            disp('Predicting LogKm half-life values (Log10 days)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogKmHL values... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.KM.Desc_i);
        
        pred = nnrpred(Xtest,train.KM.model.set.train,train.KM.model.set.y,train.KM.model.set.K,train.KM.model.set.dist_type,train.KM.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogKM_exp=NaN(size(Xtest,1),1);
        end
        res.LogKM_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.KM.model.set.train,Xtest,'auto');
        res.AD_KM=abs(AD.inorout-1)';
        res.AD_KM(round(pred.D,3)==0)=1;
        
        
        %res.AD_index=1./(1+nanmedian(pred.dc,2));
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        res.AD_index_KM=zeros(size(Xtest,1),1);
        res.Conf_index_KM=zeros(size(Xtest,1),1);
        
        KM_CAS_neighbor=cell(size(Xtest,1),5);
        KM_InChiKey_neighbor=cell(size(Xtest,1),5);
        KM_DTXSID_neighbor=cell(size(Xtest,1),5);
        KM_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogKM_Exp_neighbor=zeros(size(Xtest,1),5);
        LogKM_pred_neighbor=zeros(size(Xtest,1),5);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.KM.CAS);
                    if Li
                        res.LogKM_exp(i,1)=train.KM.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.KM.DTXSID);
                        if Li
                            res.LogKM_exp(i)=train.KM.model.set.y(Lo);
                        end
                    end
                end
            end
            KM_CAS_neighbor(i,:)=train.KM.CAS(pred.neighbors(i,:));
            KM_InChiKey_neighbor(i,:)=train.KM.InChiKey(pred.neighbors(i,:));
            KM_DTXSID_neighbor(i,:)=train.KM.DTXSID(pred.neighbors(i,:));
            KM_DSSTOXMPID_neighbor(i,:)=train.KM.DSSTOXMPID(pred.neighbors(i,:));
            LogKM_Exp_neighbor(i,:)=train.KM.model.set.y(pred.neighbors(i,:));
            LogKM_pred_neighbor(i,:)=train.KM.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogKM_Exp_neighbor(i,:),res.LogKM_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_KM(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_KM(i,1)=((1/(1+sqrt(((LogKM_Exp_neighbor(i,:)-LogKM_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_KM(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogKM_pred(i,1)=NaN;
                res.AD_KM(i)=0;
                res.AD_index_KM(i)=0;
                res.Conf_index_KM(i,1)=0;
            end
            if neighbors==1
                res.KM_CAS_neighbor(i,:)=KM_CAS_neighbor(i,:);
                res.KM_InChiKey_neighbor(i,:)=KM_InChiKey_neighbor(i,:);
                res.KM_DTXSID_neighbor(i,:)=KM_DTXSID_neighbor(i,:);
                res.KM_DSSTOXMPID_neighbor(i,:)=KM_DSSTOXMPID_neighbor(i,:);
                res.LogKM_Exp_neighbor(i,:)=LogKM_Exp_neighbor(i,:);
                res.LogKM_pred_neighbor(i,:)=LogKM_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogKM experimental= %.3f\n', res.LogKM_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogKM predicted= %.3f\n', res.LogKM_pred(i));
                if res.AD_KM(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_KM(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_KM(i));
                %CAS=strjoin(res.KM_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KM.model.set.K, res.KM_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KM.model.set.K, res.LogKM_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KM.model.set.K, res.LogKM_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogKM experimental= %.3f\n', res.LogKM_exp(i));
                end
                fprintf(output,'LogKmHL predicted= %.3f\n', res.LogKM_pred(i));
                if res.AD_KM(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_KM(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_KM(i));
                %CAS=strjoin(res.KM_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KM.model.set.K, res.KM_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KM.model.set.K, res.LogKM_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KM.model.set.K, res.LogKM_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.KM=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    %Predict KOC values
    %case {'logkoc','koc'}
    [Lia,Locb] =ismember({'koc','logkoc'},lower(prop));
    if find(Lia)
        
        
        Desc=train.KOC.Desc;
        
        
        if verbose>0
            disp('Predicting LogKoc values (Log10 L/Kg)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting LogKoc values... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        Xtest=Xin(:,train.KOC.Desc_i);
        
        pred = nnrpred(Xtest,train.KOC.model.set.train,train.KOC.model.set.y,train.KOC.model.set.K,train.KOC.model.set.dist_type,train.KOC.model.set.param.pret_type);
        pred.D=diag(pred.D);
        res.MoleculeID=MoleculeNames;
        if exp
            res.LogKoc_exp=NaN(size(Xtest,1),1);
        end
        res.LogKoc_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.KOC.model.set.train,Xtest,'auto');
        res.AD_LogKoc=abs(AD.inorout-1)';
        res.AD_LogKoc(round(pred.D,3)==0)=1;
        
        
        
        %res.AD_index=1./(1+median(pred.dc(~isnan(pred.dc)),2));
        
        
        %             res.AD_index=1./(1+median(pred.dc,2));
        %             if isnan(res.AD_index)
        %                 res.AD_index=0;
        %             end
        
        
        res.AD_index_LogKoc=zeros(size(Xtest,1),1);
        res.Conf_index_LogKoc=zeros(size(Xtest,1),1);
        
        Koc_CAS_neighbor=cell(size(Xtest,1),5);
        Koc_InChiKey_neighbor=cell(size(Xtest,1),5);
        Koc_DTXSID_neighbor=cell(size(Xtest,1),5);
        Koc_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        LogKoc_Exp_neighbor=zeros(size(Xtest,1),5);
        LogKoc_pred_neighbor=zeros(size(Xtest,1),5);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.KOC.CAS);
                    if Li
                        res.LogKoc_exp(i,1)=train.KOC.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.KOC.DTXSID);
                        if Li
                            res.LogKoc_exp(i)=train.KOC.model.set.y(Lo);
                        end
                    end
                end
            end
            Koc_CAS_neighbor(i,:)=train.KOC.CAS(pred.neighbors(i,:));
            Koc_InChiKey_neighbor(i,:)=train.KOC.InChiKey(pred.neighbors(i,:));
            Koc_DTXSID_neighbor(i,:)=train.KOC.DTXSID(pred.neighbors(i,:));
            Koc_DSSTOXMPID_neighbor(i,:)=train.KOC.DSSTOXMPID(pred.neighbors(i,:));
            LogKoc_Exp_neighbor(i,:)=train.KOC.model.set.y(pred.neighbors(i,:));
            LogKoc_pred_neighbor(i,:)=train.KOC.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.LogKoc_Exp_neighbor(i,:),res.LogKoc_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_LogKoc(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_LogKoc(i,1)=((1/(1+sqrt(((LogKoc_Exp_neighbor(i,:)-LogKoc_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_LogKoc(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.LogKoc_pred(i,1)=NaN;
                res.AD_LogKoc(i)=0;
                res.AD_index_LogKoc(i)=0;
                res.Conf_index_LogKoc(i,1)=0;
            end
            if neighbors==1
                res.Koc_CAS_neighbor(i,:)=Koc_CAS_neighbor(i,:);
                res.Koc_InChiKey_neighbor(i,:)=Koc_InChiKey_neighbor(i,:);
                res.Koc_DTXSID_neighbor(i,:)=Koc_DTXSID_neighbor(i,:);
                res.Koc_DSSTOXMPID_neighbor(i,:)=Koc_DSSTOXMPID_neighbor(i,:);
                res.LogKoc_Exp_neighbor(i,:)=LogKoc_Exp_neighbor(i,:);
                res.LogKoc_pred_neighbor(i,:)=LogKoc_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'LogKOC experimental= %.3f\n', res.LogKoc_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'LogKOC predicted= %.3f\n', res.LogKoc_pred(i));
                if res.AD_LogKoc(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_LogKoc(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_LogKoc(i));
                %CAS=strjoin(res.Koc_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KOC.model.set.K, res.Koc_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KOC.model.set.K, res.LogKoc_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KOC.model.set.K, res.LogKoc_pred_neighbor(i,1:5));
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0

                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'LogKOC experimental= %.3f\n', res.LogKoc_exp(i));
                end
                fprintf(output,'LogKOC predicted= %.3f\n', res.LogKoc_pred(i));
                if res.AD_LogKoc(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_LogKoc(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_LogKoc(i));
                %CAS=strjoin(res.Koc_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.KOC.model.set.K, res.Koc_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.KOC.model.set.K, res.LogKoc_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.KOC.model.set.K, res.LogKoc_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.KOC=res;
            clear('res');
        end
        % Clean memory
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    % ADME
    
    %Predict FUB values
    %case {'fub','fu'}
    [Lia,Locb] =ismember({'fu','fub'},lower(prop));
    if find(Lia)
        
        Desc=train.FUB.Desc;
        
        if verbose>0
            disp('Predicting FuB values (fraction)...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting FuB values... \n\n			==============================================================  \n\n');
        end
        
        %             Temp=XinCDK(:,train.FUB.cdk_in);
        %             i=1;
        %             XinCDK_FUB=zeros(size(Temp));
        %             while i<=length(train.FUB.cdk_in)
        %                 if cellfun(@ischar,table2cell(Temp(1,i)))
        %                     XinCDK_FUB(:,i)=str2double(table2cell(Temp(:,i)));
        %                 else
        %                     XinCDK_FUB(:,i)=Temp{:,i};
        %                 end
        %                 i=i+1;
        %             end
        %             clear('Temp');
        
        XinCDK_FUB=XinCDK(:,train.FUB.cdk_in);
        Xtest=[Xin(:,train.PadelVarIn(train.FUB.Padel_in)), XinCDK_FUB];
        
        Xtest=Xtest(:,train.FUB.Desc_i);
        
        pred = nnrpred(Xtest,train.FUB.model.set.train,train.FUB.model.set.y,train.FUB.model.set.K,train.FUB.model.set.dist_type,train.FUB.model.set.param.pret_type);
        pred.D=diag(pred.D);
        %pred.D=[];
        res.MoleculeID=MoleculeNames;
        if exp
            res.FUB_exp=NaN(size(Xtest,1),1);
        end
        res.FUB_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.FUB.model.set.train,Xtest,'auto');
        res.AD_FUB=abs(AD.inorout-1)';
        res.AD_FUB(round(pred.D,3)==0)=1;
        
        res.AD_index_FUB=zeros(size(Xtest,1),1);
        res.Conf_index_FUB=zeros(size(Xtest,1),1);
        
        FUB_CAS_neighbor=cell(size(Xtest,1),5);
        %FUB_InChiKey_neighbor=cell(size(Xtest,1),5);
        FUB_DTXSID_neighbor=cell(size(Xtest,1),5);
        %FUB_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        FUB_Exp_neighbor=zeros(size(Xtest,1),5);
        FUB_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.FUB.CAS);
                    if Li
                        res.FUB_exp(i,1)=train.FUB.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.FUB.DTXSID);
                        if Li
                            res.FUB_exp(i)=train.FUB.model.set.y(Lo);
                        end
                    end
                end
            end
            FUB_CAS_neighbor(i,:)=train.FUB.CAS(pred.neighbors(i,:));
            %FUB_InChiKey_neighbor(i,:)=train.FUB.InChiKey(pred.neighbors(i,:));
            FUB_DTXSID_neighbor(i,:)=train.FUB.DTXSID(pred.neighbors(i,:));
            %FUB_DSSTOXMPID_neighbor(i,:)=train.FUB.DSSTOXMPID(pred.neighbors(i,:));
            FUB_Exp_neighbor(i,:)=train.FUB.model.set.y(pred.neighbors(i,:));
            FUB_pred_neighbor(i,:)=train.FUB.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.FUB_Exp_neighbor(i,:),res.FUB_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_FUB(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_FUB(i,1)=((1/(1+sqrt(((FUB_Exp_neighbor(i,:)-FUB_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_FUB(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.FUB_pred(i,1)=NaN;
                res.AD_FUB(i)=0;
                res.AD_index_FUB(i)=0;
                res.Conf_index_FUB(i,1)=0;
            end
            
            if neighbors==1
                res.FUB_CAS_neighbor(i,:)=FUB_CAS_neighbor(i,:);
                %res.FUB_InChiKey_neighbor(i,:)=FUB_InChiKey_neighbor(i,:);
                res.FUB_DTXSID_neighbor(i,:)=FUB_DTXSID_neighbor(i,:);
                %res.FUB_DSSTOXMPID_neighbor(i,:)=FUB_DSSTOXMPID_neighbor(i,:);
                res.FUB_Exp_neighbor(i,:)=FUB_Exp_neighbor(i,:);
                res.FUB_pred_neighbor(i,:)=FUB_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'FUB experimental= %.3f\n', res.FUB_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'FUB predicted= %.3f\n', res.FUB_pred(i));
                if res.AD_FUB(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_FUB(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_FUB(i));
                %CAS=strjoin(res.FUB_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.FUB.model.set.K, res.FUB_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.FUB.model.set.K, res.FUB_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.FUB.model.set.K, res.FUB_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'FUB experimental= %.3f\n', res.FUB_exp(i));
                end
                fprintf(output,'FUB predicted= %.3f\n', res.FUB_pred(i));
                if res.AD_FUB(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_FUB(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_FUB(i));
                %CAS=strjoin(res.FUB_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.FUB.model.set.K, res.FUB_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.FUB.model.set.K, res.FUB_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.FUB.model.set.K, res.FUB_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.FUB=res;
            clear('res');
        end
        % Clean memory
        clear('XinCDK_FUB');
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    %Predict Clint values
    %case {'clint','cl'}
    [Lia,Locb] =ismember({'cl','clint'},lower(prop));
    if find(Lia)
        
        Desc=train.Clint.Desc;
        
        if verbose>0
            disp('Predicting Clint values...');
            if verbose>1
                disp(['Weighted kNN model with ', num2str(length(Desc)),'descriptors']);
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting Clint values... \n\n			==============================================================  \n\n');
        end
        
        %             Temp=XinCDK(:,train.Clint.cdk_in);
        %             i=1;
        %             XinCDK_Clint=zeros(size(Temp));
        %             while i<=length(train.Clint.cdk_in)
        %                 if cellfun(@ischar,table2cell(Temp(1,i)))
        %                     XinCDK_Clint(:,i)=str2double(table2cell(Temp(:,i)));
        %                 else
        %                     XinCDK_Clint(:,i)=Temp{:,i};
        %                 end
        %                 i=i+1;
        %             end
        %             clear('Temp');
        
        XinCDK_Clint=XinCDK(:,train.Clint.cdk_in);
        Xtest=[Xin(:,train.PadelVarIn(train.Clint.Padel_in)), XinCDK_Clint];
        
        Xtest=Xtest(:,train.Clint.Desc_i);
        
        pred = nnrpred(Xtest,train.Clint.model.set.train,train.Clint.model.set.y,train.Clint.model.set.K,train.Clint.model.set.dist_type,train.Clint.model.set.param.pret_type);
        pred.D=diag(pred.D);
        %pred.D=[];
        res.MoleculeID=MoleculeNames;
        if exp
            res.Clint_exp=NaN(size(Xtest,1),1);
        end
        res.Clint_pred(:,1)=pred.y_pred_weighted;
        AD=classical_leverage(train.Clint.model.set.train,Xtest,'auto');
        res.AD_Clint=abs(AD.inorout-1)';
        res.AD_Clint(round(pred.D,3)==0)=1;
        
        res.AD_index_Clint=zeros(size(Xtest,1),1);
        res.Conf_index_Clint=zeros(size(Xtest,1),1);
        
        Clint_CAS_neighbor=cell(size(Xtest,1),5);
        Clint_InChiKey_neighbor=cell(size(Xtest,1),5);
        Clint_DTXSID_neighbor=cell(size(Xtest,1),5);
        %Clint_DSSTOXMPID_neighbor=cell(size(Xtest,1),5);
        Clint_Exp_neighbor=zeros(size(Xtest,1),5);
        Clint_pred_neighbor=zeros(size(Xtest,1),5);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.Clint.CAS);
                    if Li
                        res.Clint_exp(i,1)=train.Clint.model.set.y(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.Clint.DTXSID);
                        if Li
                            res.Clint_exp(i)=train.Clint.model.set.y(Lo);
                        end
                    end
                end
            end
            Clint_CAS_neighbor(i,:)=train.Clint.CAS(pred.neighbors(i,:));
            Clint_InChiKey_neighbor(i,:)=train.Clint.InChiKey(pred.neighbors(i,:));
            Clint_DTXSID_neighbor(i,:)=train.Clint.DTXSID(pred.neighbors(i,:));
            %Clint_DSSTOXMPID_neighbor(i,:)=train.Clint.DSSTOXMPID(pred.neighbors(i,:));
            Clint_Exp_neighbor(i,:)=train.Clint.model.set.y(pred.neighbors(i,:));
            Clint_pred_neighbor(i,:)=train.Clint.model.yc_weighted(pred.neighbors(i,:));
            
            %                 rmse=calc_reg_param(res.Clint_Exp_neighbor(i,:),res.Clint_pred_neighbor(i,:));
            %                 res.Conf_index(i,1)=1/(1+rmse.RMSEC);
            
            res.AD_index_Clint(i,1)=1./(1+pred.dc(i,~isnan(pred.dc(i,:)))*pred.w(i,~isnan(pred.dc(i,:)))');
            
            res.Conf_index_Clint(i,1)=((1/(1+sqrt(((Clint_Exp_neighbor(i,:)-Clint_pred_neighbor(i,:)).^2)*pred.w(i,:)')))+res.AD_index_Clint(i,1))/2;
            if isempty(find(~isnan(pred.dc(i,:)), 1))
                res.Clint_pred(i,1)=NaN;
                res.AD_Clint(i)=0;
                res.AD_index_Clint(i)=0;
                res.Conf_index_Clint(i,1)=0;
            end
            if neighbors==1
                res.Clint_CAS_neighbor(i,:)=Clint_CAS_neighbor(i,:);
                res.Clint_InChiKey_neighbor(i,:)=Clint_InChiKey_neighbor(i,:);
                res.Clint_DTXSID_neighbor(i,:)=Clint_DTXSID_neighbor(i,:);
                %res.Clint_DSSTOXMPID_neighbor(i,:)=Clint_DSSTOXMPID_neighbor(i,:);
                res.Clint_Exp_neighbor(i,:)=Clint_Exp_neighbor(i,:);
                res.Clint_pred_neighbor(i,:)=Clint_pred_neighbor(i,:);
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'Clint experimental= %.3f\n', res.Clint_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'Clint predicted= %.3f\n', res.Clint_pred(i));
                if res.AD_Clint(i)==1
                    fprintf(output(Locb(find(Locb))),'AD: inside\n');
                else
                    fprintf(output(Locb(find(Locb))),'AD: outside\n');
                end
                fprintf(output(Locb(find(Locb))),'AD_index= %.2f\n', res.AD_index_Clint(i));
                fprintf(output(Locb(find(Locb))),'Conf_index= %.2f\n', res.Conf_index_Clint(i));
                %CAS=strjoin(res.Clint_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.Clint.model.set.K, res.Clint_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.Clint.model.set.K, res.Clint_Exp_neighbor(i,1:5));
                    fprintf(output(Locb(find(Locb))),'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.Clint.model.set.K, res.Clint_pred_neighbor(i,1:5));
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                %res.Xtest=Xtest;
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'Clint experimental= %.3f\n', res.Clint_exp(i));
                end
                fprintf(output,'Clint predicted= %.3f\n', res.Clint_pred(i));
                if res.AD_Clint(i)==1
                    fprintf(output,'AD: inside\n');
                else
                    fprintf(output,'AD: outside\n');
                end
                fprintf(output,'AD_index= %.2f\n', res.AD_index_Clint(i));
                fprintf(output,'Conf_index= %.2f\n', res.Conf_index_Clint(i));
                %CAS=strjoin(res.Clint_CAS_neighbor(i,1:5),',\t');
                %CAS=strrep([res.CAS_neighbors(i,1:5)],' ',', ');
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.Clint.model.set.K, res.Clint_CAS_neighbor{i,1:5});
                    fprintf(output,'Exp of the %i nearest neighbors:%15.3f,%15.3f,%15.3f,%15.3f,%15.3f\n',train.Clint.model.set.K, res.Clint_Exp_neighbor(i,1:5));
                    fprintf(output,'Pred of the %i nearest neighbors:%14.3f,%15.3f,%15.3f,%15.3f,%15.3f\n\n',train.Clint.model.set.K, res.Clint_pred_neighbor(i,1:5));
                end
                
            end
        end
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            Xtest(:,ismember(Desc,DescNames))=[];
            
            Desc(ismember(Desc,DescNames))=[];
            
            DescNames=[DescNames Desc];
            
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.Clint=res;
            clear('res');
        end
        % Clean memory
        clear('XinCDK_Clint');
        clear('Xtest');
        clear('pred');
        clear('AD');
        %end clean memory
    end
    
    % Tox properties
    
    if verbose> 0 && (tox||all)
        fprintf(1,'---------- Toxcity Endpoints ----------\n');
    end
    
    %--------------------------------------------
    
    %Predict CERAPP endpoints
    %case {'CERAPP','ER'}
    [Lia,Locb] =ismember({'cerapp','er'},lower(prop));
    if find(Lia)
        
        Desc=train.CERAPP.DescIn;
        
        if verbose>0
            disp('Predicting Estrogen Receptor Activity (CERAPP)...');
            if verbose>1
                disp('Agonist, Antagonist & Binding consensus models from the CERAPP project.');
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting CERAPP endpoints... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        
        %XlabelsCDK
        
        XinCDK_CERAPP=XinCDK(:,train.CERAPP.cdk_in);
        Xtest=[Xin(:,train.PadelVarIn(train.CERAPP.Padel_in)), XinCDK_CERAPP];
        
        XtestAG=Xtest(:,train.CERAPP.model_AG.DescAG_i);
        XtestAN=Xtest(:,train.CERAPP.model_AN.DescAN_i);
        XtestBD=Xtest(:,train.CERAPP.model_BD.DescBD_i);
        
        predAG = knnpred2(XtestAG,train.CERAPP.model_AG.set.train,train.CERAPP.model_AG.set.class,train.CERAPP.model_AG.set.K,train.CERAPP.model_AG.set.dist_type,train.CERAPP.model_AG.set.param.pret_type);
        %predAG.D=diag(predAG.D);
        predAG.D=[];
        predAN = knnpred2(XtestAN,train.CERAPP.model_AN.set.train,train.CERAPP.model_AN.set.class,train.CERAPP.model_AN.set.K,train.CERAPP.model_AN.set.dist_type,train.CERAPP.model_AN.set.param.pret_type);
        %predAN.D=diag(predAN.D);
        predAN.D=[];
        predBD = knnpred2(XtestBD,train.CERAPP.model_BD.set.train,train.CERAPP.model_BD.set.class,train.CERAPP.model_BD.set.K,train.CERAPP.model_BD.set.dist_type,train.CERAPP.model_BD.set.param.pret_type);
        %predBD.D=diag(predBD.D);
        predBD.D=[];
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.CERAPP_Ago_exp=cell(size(Xtest,1),1);
        end
        res.CERAPP_Ago_pred(:,1)=predAG.class_pred-1;
        AD=classical_leverage(train.CERAPP.model_AG.set.train,XtestAG,'auto');
        res.AD_CERAPP_Ago=abs(AD.inorout-1)';
        res.AD_index_CERAPP_Ago=1-test_pretreatment(predAG.dc(:,1),train.CERAPP.model_AG.set.dc_param);
        res.AD_index_CERAPP_Ago(find(res.AD_index_CERAPP_Ago<0),1)=1./(1+predAG.dc(find(res.AD_index_CERAPP_Ago<0),1));
        res.AD_CERAPP_Ago(find(isnan(predAG.dc(:,1))))=0;
        res.AD_index_CERAPP_Ago(find(isnan(predAG.dc(:,1))))=0;
        res.AD_CERAPP_Ago(find(res.AD_index_CERAPP_Ago>0.5))=1;
        res.Conf_index_CERAPP_Ago=zeros(size(XtestAG,1),1);
        if exp
            res.CERAPP_Anta_exp=cell(size(Xtest,1),1);
        end
        res.CERAPP_Anta_pred(:,1)=predAN.class_pred-1;
        AD=classical_leverage(train.CERAPP.model_AN.set.train,XtestAN,'auto');
        res.AD_CERAPP_Anta=abs(AD.inorout-1)';
        res.AD_index_CERAPP_Anta=1-test_pretreatment(predAN.dc(:,1),train.CERAPP.model_AN.set.dc_param);
        res.AD_index_CERAPP_Anta(find(res.AD_index_CERAPP_Anta<0),1)=1./(1+predAN.dc(find(res.AD_index_CERAPP_Anta<0),1));
        res.AD_CERAPP_Anta(find(isnan(predAN.dc(:,1))))=0;
        res.AD_index_CERAPP_Anta(find(isnan(predAN.dc(:,1))))=0;
        res.AD_CERAPP_Anta(find(res.AD_index_CERAPP_Anta>0.5))=1;
        res.Conf_index_CERAPP_Anta=zeros(size(XtestAN,1),1);
        if exp
            res.CERAPP_Bind_exp=cell(size(Xtest,1),1);
        end
        res.CERAPP_Bind_pred(:,1)=predBD.class_pred-1;
        AD=classical_leverage(train.CERAPP.model_BD.set.train,XtestBD,'auto');
        res.AD_CERAPP_Bind=abs(AD.inorout-1)';
        res.AD_index_CERAPP_Bind=1-test_pretreatment(predBD.dc(:,1),train.CERAPP.model_BD.set.dc_param);
        res.AD_index_CERAPP_Bind(find(res.AD_index_CERAPP_Bind<0),1)=1./(1+predBD.dc(find(res.AD_index_CERAPP_Bind<0),1));
        res.AD_CERAPP_Bind(find(isnan(predBD.dc(:,1))))=0;
        res.AD_index_CERAPP_Bind(find(isnan(predBD.dc(:,1))))=0;
        res.AD_CERAPP_Bind(find(res.AD_index_CERAPP_Bind>0.5))=1;
        res.Conf_index_CERAPP_Bind=zeros(size(XtestBD,1),1);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.CERAPP.model_AG.CAS);
                    if Li
                        res.CERAPP_Ago_exp(i,1)=train.CERAPP.model_AG.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CERAPP.model_AG.DTXSID);
                        if Li
                            res.CERAPP_Ago_exp(i)=train.CERAPP.model_AG.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CERAPP.model_AN.CAS);
                    if Li
                        res.CERAPP_Anta_exp(i,1)=train.CERAPP.model_AN.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CERAPP.model_AN.DTXSID);
                        if Li
                            res.CERAPP_Anta_exp(i)=train.CERAPP.model_AN.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CERAPP.model_BD.CAS);
                    if Li
                        res.CERAPP_Bind_exp(i,1)=train.CERAPP.model_BD.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CERAPP.model_BD.DTXSID);
                        if Li
                            res.CERAPP_Bind_exp(i)=train.CERAPP.model_BD.set.class_Exp(Lo);
                        end
                    end
                end
            end
            res.Conf_index_CERAPP_Ago(i,1)=train.CERAPP.model_AG.conc_AG(predAG.neighbors(i,:),1)'*predAG.w(i,:)';
            res.Conf_index_CERAPP_Anta(i,1)=train.CERAPP.model_AN.conc_AN(predAN.neighbors(i,:),1)'*predAN.w(i,:)';
            res.Conf_index_CERAPP_Bind(i,1)=train.CERAPP.model_BD.conc_BD(predBD.neighbors(i,:),1)'*predBD.w(i,:)';
            
            if neighbors==1
                res.CERAPP_Ago_CAS_neighbor(i,:)=train.CERAPP.model_AG.CAS(predAG.neighbors(i,:));
                res.CERAPP_Ago_InChiKey_neighbor(i,:)=train.CERAPP.model_AG.InChiKey(predAG.neighbors(i,:));
                res.CERAPP_Ago_DTXSID_neighbor(i,:)=train.CERAPP.model_AG.DTXSID(predAG.neighbors(i,:));
                %res.CERAPP_Ago_DSSTOXMPID_neighbor(i,:)=train.CERAPP.model_AG.DSSTOXMPID(pred.neighbors(i,:));
                res.CERAPP_Ago_Exp_neighbor(i,:)=train.CERAPP.model_AG.set.class_Exp(predAG.neighbors(i,:));
                res.CERAPP_Ago_pred_neighbor(i,:)=train.CERAPP.model_AG.set.class_S(predAG.neighbors(i,:));
                
                res.CERAPP_Anta_CAS_neighbor(i,:)=train.CERAPP.model_AN.CAS(predAN.neighbors(i,:));
                res.CERAPP_Anta_InChiKey_neighbor(i,:)=train.CERAPP.model_AN.InChiKey(predAN.neighbors(i,:));
                res.CERAPP_Anta_DTXSID_neighbor(i,:)=train.CERAPP.model_AN.DTXSID(predAN.neighbors(i,:));
                %res.CERAPP_Anta_DSSTOXMPID_neighbor(i,:)=train.CERAPP.model_AN.DSSTOXMPID(pred.neighbors(i,:));
                res.CERAPP_Anta_Exp_neighbor(i,:)=train.CERAPP.model_AN.set.class_Exp(predAN.neighbors(i,:));
                res.CERAPP_Anta_pred_neighbor(i,:)=train.CERAPP.model_AN.set.class_S(predAN.neighbors(i,:));
                
                res.CERAPP_Bind_CAS_neighbor(i,:)=train.CERAPP.model_BD.CAS(predBD.neighbors(i,:));
                res.CERAPP_Bind_InChiKey_neighbor(i,:)=train.CERAPP.model_BD.InChiKey(predBD.neighbors(i,:));
                res.CERAPP_Bind_DTXSID_neighbor(i,:)=train.CERAPP.model_BD.DTXSID(predBD.neighbors(i,:));
                %res.CERAPP_Bind_DSSTOXMPID_neighbor(i,:)=train.CERAPP.model_BD.DSSTOXMPID(predBD.neighbors(i,:));
                res.CERAPP_Bind_Exp_neighbor(i,:)=train.CERAPP.model_BD.set.class_Exp(predBD.neighbors(i,:));
                res.CERAPP_Bind_pred_neighbor(i,:)=train.CERAPP.model_BD.set.class_S(predBD.neighbors(i,:));
                
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'AG experimental= %s, AN experimental= %s, BD experimental= %s\n', res.CERAPP_Ago_exp{i},res.CERAPP_Anta_exp{i},res.CERAPP_Bind_exp{i});
                end
                fprintf(output(Locb(find(Locb))),'AG predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Ago_pred(i),res.AD_CERAPP_Ago(i),res.AD_index_CERAPP_Ago(i),res.Conf_index_CERAPP_Ago(i));
                fprintf(output(Locb(find(Locb))),'AN predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Anta_pred(i),res.AD_CERAPP_Anta(i),res.AD_index_CERAPP_Anta(i),res.Conf_index_CERAPP_Anta(i));
                fprintf(output(Locb(find(Locb))),'BD predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Bind_pred(i),res.AD_CERAPP_Bind(i),res.AD_index_CERAPP_Bind(i),res.Conf_index_CERAPP_Bind(i));
                
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_AG.set.K, res.CERAPP_Ago_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_AN.set.K, res.CERAPP_Anta_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_BD.set.K, res.CERAPP_Bind_CAS_neighbor{i,1:5});
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'AG experimental= %s, AN experimental= %s, BD experimental= %s\n', res.CERAPP_Ago_exp{i},res.CERAPP_Anta_exp{i},res.CERAPP_Bind_exp{i});
                end
                fprintf(output,'AG predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Ago_pred(i),res.AD_CERAPP_Ago(i),res.AD_index_CERAPP_Ago(i),res.Conf_index_CERAPP_Ago(i));
                fprintf(output,'AN predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Anta_pred(i),res.AD_CERAPP_Anta(i),res.AD_index_CERAPP_Anta(i),res.Conf_index_CERAPP_Anta(i));
                fprintf(output,'BD predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CERAPP_Bind_pred(i),res.AD_CERAPP_Bind(i),res.AD_index_CERAPP_Bind(i),res.Conf_index_CERAPP_Bind(i));
                
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_AG.set.K, res.CERAPP_Ago_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_AN.set.K, res.CERAPP_Anta_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CERAPP.model_BD.set.K, res.CERAPP_Bind_CAS_neighbor{i,1:5});
                end
                
            end
        end
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                %Xtest=[XtestAG; XtestAN; XtestBD; XtestGHS; XtestLD50];
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            %Xtest=[XtestAG; XtestAN; XtestBD; XtestGHS; XtestLD50];
            Xtest(:,ismember(Desc,DescNames))=[];
            Desc(ismember(Desc,DescNames))=[];
            DescNames=[DescNames Desc];
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.CERAPP=res;
            clear('res');
            
        end
        % Clean memory
        clear('XinCDK_CERAPP');
        clear('Xtest');
        clear('XtestAG');
        clear('XtestAN');
        clear('XtestBD');
        clear('predAG');
        clear('predAN');
        clear('predBD')
        clear('AD');
        %end clean memory
    end
    
    %--------------------------------------------
    
    %Predict CoMPARA endpoints
    %case {'CoMPARA','AR'}
    [Lia,Locb] =ismember({'compara','ar'},lower(prop));
    if find(Lia)
        
        Desc=train.CoMPARA.DescIn;
        
        if verbose>0
            disp('Predicting Androgen Receptor Activity (CoMPARA)...');
            if verbose>1
                disp('Agonist, Antagonist & Binding consensus models from the CATMoS project.');
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting CoMPARA endpoints... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        
        %XlabelsCDK
        XinCDK_CoMPARA=XinCDK(:,train.CoMPARA.cdk_in);
        Xtest=[Xin(:,train.PadelVarIn(train.CoMPARA.Padel_in)), XinCDK_CoMPARA];
        
        XtestAG=Xtest(:,train.CoMPARA.model_AG.DescAG_i);
        XtestAN=Xtest(:,train.CoMPARA.model_AN.DescAN_i);
        XtestBD=Xtest(:,train.CoMPARA.model_BD.DescBD_i);
        
        predAG = knnpred2(XtestAG,train.CoMPARA.model_AG.set.train,train.CoMPARA.model_AG.set.class,train.CoMPARA.model_AG.set.K,train.CoMPARA.model_AG.set.dist_type,train.CoMPARA.model_AG.set.param.pret_type);
        %predAG.D=diag(predAG.D);
        predAG.D=[];
        predAN = knnpred2(XtestAN,train.CoMPARA.model_AN.set.train,train.CoMPARA.model_AN.set.class,train.CoMPARA.model_AN.set.K,train.CoMPARA.model_AN.set.dist_type,train.CoMPARA.model_AN.set.param.pret_type);
        %predAN.D=diag(predAN.D);
        predAN.D=[];
        predBD = knnpred2(XtestBD,train.CoMPARA.model_BD.set.train,train.CoMPARA.model_BD.set.class,train.CoMPARA.model_BD.set.K,train.CoMPARA.model_BD.set.dist_type,train.CoMPARA.model_BD.set.param.pret_type);
        %predBD.D=diag(predBD.D);
        predBD.D=[];
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.CoMPARA_Ago_exp=cell(size(Xtest,1),1);
        end
        res.CoMPARA_Ago_pred(:,1)=predAG.class_pred-1;
        AD=classical_leverage(train.CoMPARA.model_AG.set.train,XtestAG,'auto');
        res.AD_CoMPARA_Ago=abs(AD.inorout-1)';
        res.AD_index_CoMPARA_Ago=1-test_pretreatment(predAG.dc(:,1),train.CoMPARA.model_AG.set.dc_param);
        res.AD_index_CoMPARA_Ago(find(res.AD_index_CoMPARA_Ago<0),1)=1./(1+predAG.dc(find(res.AD_index_CoMPARA_Ago<0),1));
        res.AD_CoMPARA_Ago(find(isnan(predAG.dc(:,1))))=0;
        res.AD_index_CoMPARA_Ago(find(isnan(predAG.dc(:,1))))=0;
        res.AD_CoMPARA_Ago(find(res.AD_index_CoMPARA_Ago>0.5))=1;
        res.Conf_index_CoMPARA_Ago=zeros(size(XtestAG,1),1);
        if exp
            res.CoMPARA_Anta_exp=cell(size(Xtest,1),1);
        end
        res.CoMPARA_Anta_pred(:,1)=predAN.class_pred-1;
        AD=classical_leverage(train.CoMPARA.model_AN.set.train,XtestAN,'auto');
        res.AD_CoMPARA_Anta=abs(AD.inorout-1)';
        res.AD_index_CoMPARA_Anta=1-test_pretreatment(predAN.dc(:,1),train.CoMPARA.model_AN.set.dc_param);
        res.AD_index_CoMPARA_Anta(find(res.AD_index_CoMPARA_Anta<0),1)=1./(1+predAN.dc(find(res.AD_index_CoMPARA_Anta<0),1));
        res.AD_CoMPARA_Anta(find(isnan(predAN.dc(:,1))))=0;
        res.AD_index_CoMPARA_Anta(find(isnan(predAN.dc(:,1))))=0;
        res.AD_CoMPARA_Anta(find(res.AD_index_CoMPARA_Anta>0.5))=1;
        res.Conf_index_CoMPARA_Anta=zeros(size(XtestAN,1),1);
        if exp
            res.CoMPARA_Bind_exp=cell(size(Xtest,1),1);
        end
        res.CoMPARA_Bind_pred(:,1)=predBD.class_pred-1;
        AD=classical_leverage(train.CoMPARA.model_BD.set.train,XtestBD,'auto');
        res.AD_CoMPARA_Bind=abs(AD.inorout-1)';
        res.AD_index_CoMPARA_Bind=1-test_pretreatment(predBD.dc(:,1),train.CoMPARA.model_BD.set.dc_param);
        res.AD_index_CoMPARA_Bind(find(res.AD_index_CoMPARA_Bind<0),1)=1./(1+predBD.dc(find(res.AD_index_CoMPARA_Bind<0),1));
        res.AD_CoMPARA_Bind(find(isnan(predBD.dc(:,1))))=0;
        res.AD_index_CoMPARA_Bind(find(isnan(predBD.dc(:,1))))=0;
        res.AD_CoMPARA_Bind(find(res.AD_index_CoMPARA_Bind>0.5))=1;
        res.Conf_index_CoMPARA_Bind=zeros(size(XtestBD,1),1);
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.CoMPARA.model_AG.CAS);
                    if Li
                        res.CoMPARA_Ago_exp(i)=train.CoMPARA.model_AG.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CoMPARA.model_AG.DTXSID);
                        if Li
                            res.CoMPARA_Ago_exp(i)=train.CoMPARA.model_AG.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CoMPARA.model_AN.CAS);
                    if Li
                        res.CoMPARA_Anta_exp(i)=train.CoMPARA.model_AN.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CoMPARA.model_AN.DTXSID);
                        if Li
                            res.CoMPARA_Anta_exp(i)=train.CoMPARA.model_AN.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CoMPARA.model_BD.CAS);
                    if Li
                        res.CoMPARA_Bind_exp(i)=train.CoMPARA.model_BD.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CoMPARA.model_BD.DTXSID);
                        if Li
                            res.CoMPARA_Bind_exp(i)=train.CoMPARA.model_BD.set.class_Exp(Lo);
                        end
                    end
                end
            end
            res.Conf_index_CoMPARA_Ago(i,1)=train.CoMPARA.model_AG.conc_AG(predAG.neighbors(i,:),1)'*predAG.w(i,:)';
            res.Conf_index_CoMPARA_Anta(i,1)=train.CoMPARA.model_AN.conc_AN(predAN.neighbors(i,:),1)'*predAN.w(i,:)';
            res.Conf_index_CoMPARA_Bind(i,1)=train.CoMPARA.model_BD.conc_BD(predBD.neighbors(i,:),1)'*predBD.w(i,:)';
            
            if neighbors==1
                res.CoMPARA_Ago_CAS_neighbor(i,:)=train.CoMPARA.model_AG.CAS(predAG.neighbors(i,:));
                res.CoMPARA_Ago_InChiKey_neighbor(i,:)=train.CoMPARA.model_AG.InChiKey(predAG.neighbors(i,:));
                res.CoMPARA_Ago_DTXSID_neighbor(i,:)=train.CoMPARA.model_AG.DTXSID(predAG.neighbors(i,:));
                %res.CoMPARA_Ago_DSSTOXMPID_neighbor(i,:)=train.CoMPARA.model_AG.DSSTOXMPID(pred.neighbors(i,:));
                res.CoMPARA_Ago_Exp_neighbor(i,:)=train.CoMPARA.model_AG.set.class_Exp(predAG.neighbors(i,:));
                res.CoMPARA_Ago_pred_neighbor(i,:)=train.CoMPARA.model_AG.set.class_S(predAG.neighbors(i,:));
                
                res.CoMPARA_Anta_CAS_neighbor(i,:)=train.CoMPARA.model_AN.CAS(predAN.neighbors(i,:));
                res.CoMPARA_Anta_InChiKey_neighbor(i,:)=train.CoMPARA.model_AN.InChiKey(predAN.neighbors(i,:));
                res.CoMPARA_Anta_DTXSID_neighbor(i,:)=train.CoMPARA.model_AN.DTXSID(predAN.neighbors(i,:));
                %res.CoMPARA_Anta_DSSTOXMPID_neighbor(i,:)=train.CoMPARA.model_AN.DSSTOXMPID(pred.neighbors(i,:));
                res.CoMPARA_Anta_Exp_neighbor(i,:)=train.CoMPARA.model_AN.set.class_Exp(predAN.neighbors(i,:));
                res.CoMPARA_Anta_pred_neighbor(i,:)=train.CoMPARA.model_AN.set.class_S(predAN.neighbors(i,:));
                
                res.CoMPARA_Bind_CAS_neighbor(i,:)=train.CoMPARA.model_BD.CAS(predBD.neighbors(i,:));
                res.CoMPARA_Bind_InChiKey_neighbor(i,:)=train.CoMPARA.model_BD.InChiKey(predBD.neighbors(i,:));
                res.CoMPARA_Bind_DTXSID_neighbor(i,:)=train.CoMPARA.model_BD.DTXSID(predBD.neighbors(i,:));
                %res.CoMPARA_Bind_DSSTOXMPID_neighbor(i,:)=train.CoMPARA.model_BD.DSSTOXMPID(predBD.neighbors(i,:));
                res.CoMPARA_Bind_Exp_neighbor(i,:)=train.CoMPARA.model_BD.set.class_Exp(predBD.neighbors(i,:));
                res.CoMPARA_Bind_pred_neighbor(i,:)=train.CoMPARA.model_BD.set.class_S(predBD.neighbors(i,:));
                
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'AG experimental= %s, AN experimental= %s, BD experimental= %s\n', res.CoMPARA_Ago_exp{i},res.CoMPARA_Anta_exp{i},res.CoMPARA_Bind_exp{i});
                end
                fprintf(output(Locb(find(Locb))),'AG predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Ago_pred(i),res.AD_CoMPARA_Ago(i),res.AD_index_CoMPARA_Ago(i),res.Conf_index_CoMPARA_Ago(i));
                fprintf(output(Locb(find(Locb))),'AN predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Anta_pred(i),res.AD_CoMPARA_Anta(i),res.AD_index_CoMPARA_Anta(i),res.Conf_index_CoMPARA_Anta(i));
                fprintf(output(Locb(find(Locb))),'BD category predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Bind_pred(i),res.AD_CoMPARA_Bind(i),res.AD_index_CoMPARA_Bind(i),res.Conf_index_CoMPARA_Bind(i));
                
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_AG.set.K, res.CoMPARA_Ago_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_AN.set.K, res.CoMPARA_Anta_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_BD.set.K, res.CoMPARA_Bind_CAS_neighbor{i,1:5});
                end
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'AG experimental= %s, AN experimental= %s, BD experimental= %s\n', res.CoMPARA_Ago_exp{i},res.CoMPARA_Anta_exp{i},res.CoMPARA_Bind_exp{i});
                end
                fprintf(output,'AG predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Ago_pred(i),res.AD_CoMPARA_Ago(i),res.AD_index_CoMPARA_Ago(i),res.Conf_index_CoMPARA_Ago(i));
                fprintf(output,'AN predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Anta_pred(i),res.AD_CoMPARA_Anta(i),res.AD_index_CoMPARA_Anta(i),res.Conf_index_CoMPARA_Anta(i));
                fprintf(output,'BD category predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CoMPARA_Bind_pred(i),res.AD_CoMPARA_Bind(i),res.AD_index_CoMPARA_Bind(i),res.Conf_index_CoMPARA_Bind(i));
                
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_AG.set.K, res.CoMPARA_Ago_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_AN.set.K, res.CoMPARA_Anta_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CoMPARA.model_BD.set.K, res.CoMPARA_Bind_CAS_neighbor{i,1:5});
                end
                
            end
        end
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                %Xtest=[XtestAG; XtestAN; XtestBD; XtestGHS; XtestLD50];
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            %Xtest=[XtestAG; XtestAN; XtestBD; XtestGHS; XtestLD50];
            Xtest(:,ismember(Desc,DescNames))=[];
            Desc(ismember(Desc,DescNames))=[];
            DescNames=[DescNames Desc];
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.CoMPARA=res;
            clear('res');
            
        end
        % Clean memory
        clear('XinCDK_CoMPARA');
        clear('Xtest');
        clear('XtestAG');
        clear('XtestAN');
        clear('XtestBD');
        clear('predAG');
        clear('predAN');
        clear('predBD')
        clear('AD');
        %end clean memory
    end
    
    %--------------------------------------------
    %Predict CATMoS endpoints
    %case {'CATMoS','AcuteTox'}
    [Lia,Locb] =ismember({'catmos','acutetox'},lower(prop));
    if find(Lia)
        
        Desc=train.CATMoS.DescIn;
        
        if verbose>0
            disp('Predicting Acute Oral Tox. endpoints (CATMoS)...');
            if verbose>1
                disp('VT, NT, EPA, GHS & LD50 consensus models from the CATMoS project.');
            end
            
        end
        
        if strcmpi(ext,'.txt') && sep==0
            fprintf(output,'\n\n\t\t\t\t\t Predicting CATMoS endpoints... \n\n			==============================================================  \n\n');
        end
        
        %             Xtest=zeros(size(Xin,1),length(Desc));
        %
        %             for i=1:length(Desc)
        %                 for l=1:length(Xin(1,:))
        %                     if strcmp(Desc(i),Xlabels(l))
        %                         Xtest(:,i)=Xin(:,l);
        %                         break;
        %                     end
        %                 end
        %             end
        
        %XlabelsCDK
        XinCDK_CATMoS=XinCDK(:,train.CATMoS.cdk_in);
        Xtest=[Xin(:,train.PadelVarIn(train.CATMoS.Padel_in)), XinCDK_CATMoS];
        
        XtestVT=Xtest(:,train.CATMoS.model_VT.DescVT_i);
        XtestNT=Xtest(:,train.CATMoS.model_NT.DescNT_i);
        XtestEPA=Xtest(:,train.CATMoS.model_EPA.DescEPA_i);
        XtestGHS=Xtest(:,train.CATMoS.model_GHS.DescGHS_i);
        XtestLD50=Xtest(:,train.CATMoS.model_LD50.DescLD50_i);
        
        predVT = knnpred2(XtestVT,train.CATMoS.model_VT.set.train,train.CATMoS.model_VT.set.class,train.CATMoS.model_VT.set.K,train.CATMoS.model_VT.set.dist_type,train.CATMoS.model_VT.set.param.pret_type);
        %predVT.D=diag(predVT.D);
        predVT.D=[];
        predNT = knnpred2(XtestNT,train.CATMoS.model_NT.set.train,train.CATMoS.model_NT.set.class,train.CATMoS.model_NT.set.K,train.CATMoS.model_NT.set.dist_type,train.CATMoS.model_NT.set.param.pret_type);
        %predNT.D=diag(predNT.D);
        predNT.D=[];
        predEPA = knnpred2(XtestEPA,train.CATMoS.model_EPA.set.train,train.CATMoS.model_EPA.set.class,train.CATMoS.model_EPA.set.K,train.CATMoS.model_EPA.set.dist_type,train.CATMoS.model_EPA.set.param.pret_type);
        %predEPA.D=diag(predEPA.D);
        predEPA.D=[];
        predGHS = knnpred2(XtestGHS,train.CATMoS.model_GHS.set.train,train.CATMoS.model_GHS.set.class,train.CATMoS.model_GHS.set.K,train.CATMoS.model_GHS.set.dist_type,train.CATMoS.model_GHS.set.param.pret_type);
        %predGHS.D=diag(predGHS.D);
        predGHS.D=[];
        predLD50 = nnrpred2(XtestLD50,train.CATMoS.model_LD50.set.train,train.CATMoS.model_LD50.set.y,train.CATMoS.model_LD50.set.K,train.CATMoS.model_LD50.set.dist_type,train.CATMoS.model_LD50.set.param.pret_type);
        %predLD50.D=diag(predLD50.D);
        predLD50.D=[];
        
        res.MoleculeID=MoleculeNames;
        if exp
            res.CATMoS_VT_exp=NaN(size(Xtest,1),1);
        end
        res.CATMoS_VT_pred(:,1)=predVT.class_pred-1;
        AD=classical_leverage(train.CATMoS.model_VT.set.train,XtestVT,'auto');
        res.AD_VT=abs(AD.inorout-1)';
        res.AD_index_VT=1-test_pretreatment(predVT.dc(:,1),train.CATMoS.model_VT.set.dc_param);
        res.AD_index_VT(find(res.AD_index_VT<0),1)=1./(1+predVT.dc(find(res.AD_index_VT<0),1));
        res.AD_VT(find(isnan(predVT.dc(:,1))))=0;
        res.AD_index_VT(find(isnan(predVT.dc(:,1))))=0;
        res.AD_VT(find(res.AD_index_VT>0.5))=1;
        res.Conf_index_VT=zeros(size(XtestVT,1),1);
        if exp
            res.CATMoS_NT_exp=NaN(size(Xtest,1),1);
        end
        res.CATMoS_NT_pred(:,1)=predNT.class_pred-1;
        AD=classical_leverage(train.CATMoS.model_NT.set.train,XtestNT,'auto');
        res.AD_NT=abs(AD.inorout-1)';
        res.AD_index_NT=1-test_pretreatment(predNT.dc(:,1),train.CATMoS.model_NT.set.dc_param);
        res.AD_index_NT(find(res.AD_index_NT<0),1)=1./(1+predNT.dc(find(res.AD_index_NT<0),1));
        res.AD_NT(find(isnan(predNT.dc(:,1))))=0;
        res.AD_index_NT(find(isnan(predNT.dc(:,1))))=0;
        res.AD_NT(find(res.AD_index_NT>0.5))=1;
        res.Conf_index_NT=zeros(size(XtestNT,1),1);
        if exp
            res.CATMoS_EPA_exp=NaN(size(Xtest,1),1);
        end
        res.CATMoS_EPA_pred(:,1)=predEPA.class_pred;
        AD=classical_leverage(train.CATMoS.model_EPA.set.train,XtestEPA,'auto');
        res.AD_EPA=abs(AD.inorout-1)';
        res.AD_index_EPA=1-test_pretreatment(predEPA.dc(:,1),train.CATMoS.model_EPA.set.dc_param);
        res.AD_index_EPA(find(res.AD_index_EPA<0),1)=1./(1+predEPA.dc(find(res.AD_index_EPA<0),1));
        res.AD_EPA(find(isnan(predEPA.dc(:,1))))=0;
        res.AD_index_EPA(find(isnan(predEPA.dc(:,1))))=0;
        res.AD_EPA(find(res.AD_index_EPA>0.5))=1;
        res.Conf_index_EPA=zeros(size(XtestEPA,1),1);
        if exp
            res.CATMoS_GHS_exp=NaN(size(Xtest,1),1);
        end
        res.CATMoS_GHS_pred(:,1)=predGHS.class_pred;
        AD=classical_leverage(train.CATMoS.model_GHS.set.train,XtestGHS,'auto');
        res.AD_GHS=abs(AD.inorout-1)';
        res.AD_index_GHS=1-test_pretreatment(predGHS.dc(:,1),train.CATMoS.model_GHS.set.dc_param);
        res.AD_index_GHS(find(res.AD_index_GHS<0),1)=1./(1+predGHS.dc(find(res.AD_index_GHS<0),1));
        res.AD_GHS(find(isnan(predGHS.dc(:,1))))=0;
        res.AD_index_GHS(find(isnan(predGHS.dc(:,1))))=0;
        res.AD_GHS(find(res.AD_index_GHS>0.5))=1;
        res.Conf_index_GHS=zeros(size(XtestGHS,1),1);
        if exp
            res.CATMoS_LD50_exp=NaN(size(Xtest,1),1);
        end
        res.CATMoS_LD50_pred(:,1)=predLD50.y_pred_weighted;
        AD=classical_leverage(train.CATMoS.model_LD50.set.train,XtestLD50,'auto');
        res.AD_LD50=abs(AD.inorout-1)';
        res.AD_index_LD50=1-test_pretreatment(predLD50.dc(:,1),train.CATMoS.model_LD50.set.dc_param);
        res.AD_index_LD50(find(res.AD_index_LD50<0),1)=1./(1+predLD50.dc(find(res.AD_index_LD50<0),1));
        res.AD_LD50(find(isnan(predLD50.dc(:,1))))=0;
        res.AD_index_LD50(find(isnan(predLD50.dc(:,1))))=0;
        res.AD_LD50(find(res.AD_index_LD50>0.5))=1;
        res.Conf_index_LD50=zeros(size(XtestLD50,1),1);
        
        
        for i=1:size(Xtest(:,1))
            if exp
                if ~contains(MoleculeNames(i),'AUTOGEN_')
                    [Li,Lo] = ismember(MoleculeNames(i),train.CATMoS.model_VT.CAS);
                    if Li
                        res.CATMoS_VT_exp(i,1)=train.CATMoS.model_VT.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CATMoS.model_VT.DTXSID);
                        if Li
                            res.CATMoS_VT_exp(i)=train.CATMoS.model_VT.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CATMoS.model_NT.CAS);
                    if Li
                        res.CATMoS_NT_exp(i,1)=train.CATMoS.model_NT.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CATMoS.model_NT.DTXSID);
                        if Li
                            res.CATMoS_NT_exp(i)=train.CATMoS.model_NT.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CATMoS.model_EPA.CAS);
                    if Li
                        res.CATMoS_EPA_exp(i,1)=train.CATMoS.model_EPA.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CATMoS.model_EPA.DTXSID);
                        if Li
                            res.CATMoS_EPA_exp(i)=train.CATMoS.model_EPA.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CATMoS.model_GHS.CAS);
                    if Li
                        res.CATMoS_GHS_exp(i,1)=train.CATMoS.model_GHS.set.class_Exp(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CATMoS.model_GHS.DTXSID);
                        if Li
                            res.CATMoS_GHS_exp(i)=train.CATMoS.model_GHS.set.class_Exp(Lo);
                        end
                    end
                    [Li,Lo] = ismember(MoleculeNames(i),train.CATMoS.model_LD50.CAS);
                    if Li
                        res.CATMoS_LD50_exp(i,1)=train.CATMoS.model_LD50.set.y_Exp_n(Lo);
                    else
                        [Li,Lo] = ismember(MoleculeNames{i},train.CATMoS.model_LD50.DTXSID);
                        if Li
                            res.CATMoS_LD50_exp(i)=train.CATMoS.model_LD50.set.y_Exp_n(Lo);
                        end
                    end
                end
            end
            %res.AD_index_VT(i,1)=1./(1+predVT.dc(i,1)*predVT.w(i,~isnan(predVT.dc(i,1)))');
            
            res.Conf_index_VT(i,1)=train.CATMoS.model_VT.conc_VT(predVT.neighbors(i,:),1)'*predVT.w(i,:)';
            res.Conf_index_NT(i,1)=train.CATMoS.model_NT.conc_NT(predNT.neighbors(i,:),1)'*predNT.w(i,:)';
            res.Conf_index_EPA(i,1)=train.CATMoS.model_EPA.conc_EPA(predEPA.neighbors(i,:),1)'*predEPA.w(i,:)';
            res.Conf_index_GHS(i,1)=train.CATMoS.model_GHS.conc_GHS(predGHS.neighbors(i,:),1)'*predGHS.w(i,:)';
            res.Conf_index_LD50(i,1)=train.CATMoS.model_LD50.conc_LD50(predLD50.neighbors(i,:),1)'*predLD50.w(i,:)';
            
            if neighbors==1
                res.VT_CATMoS_ID_neighbor(i,:)=train.CATMoS.model_VT.ChemID(predVT.neighbors(i,:));
                res.VT_CAS_neighbor(i,:)=train.CATMoS.model_VT.CAS(predVT.neighbors(i,:));
                res.VT_InChiKey_neighbor(i,:)=train.CATMoS.model_VT.InChiKey(predVT.neighbors(i,:));
                res.VT_DTXSID_neighbor(i,:)=train.CATMoS.model_VT.DTXSID(predVT.neighbors(i,:));
                %res.VT_DSSTOXMPID_neighbor(i,:)=train.CATMoS.model_VT.DSSTOXMPID(pred.neighbors(i,:));
                res.VT_Exp_neighbor(i,:)=train.CATMoS.model_VT.set.class_Exp(predVT.neighbors(i,:));
                res.VT_pred_neighbor(i,:)=train.CATMoS.model_VT.set.class(predVT.neighbors(i,:))-1;
                
                res.NT_CATMoS_ID_neighbor(i,:)=train.CATMoS.model_NT.ChemID(predNT.neighbors(i,:));
                res.NT_CAS_neighbor(i,:)=train.CATMoS.model_NT.CAS(predNT.neighbors(i,:));
                res.NT_InChiKey_neighbor(i,:)=train.CATMoS.model_NT.InChiKey(predNT.neighbors(i,:));
                res.NT_DTXSID_neighbor(i,:)=train.CATMoS.model_NT.DTXSID(predNT.neighbors(i,:));
                %res.NT_DSSTOXMPID_neighbor(i,:)=train.CATMoS.model_NT.DSSTOXMPID(pred.neighbors(i,:));
                res.NT_Exp_neighbor(i,:)=train.CATMoS.model_NT.set.class_Exp(predNT.neighbors(i,:));
                res.NT_pred_neighbor(i,:)=train.CATMoS.model_NT.set.class(predNT.neighbors(i,:))-1;
                
                res.EPA_CATMoS_ID_neighbor(i,:)=train.CATMoS.model_EPA.ChemID(predEPA.neighbors(i,:));
                res.EPA_CAS_neighbor(i,:)=train.CATMoS.model_EPA.CAS(predEPA.neighbors(i,:));
                res.EPA_InChiKey_neighbor(i,:)=train.CATMoS.model_EPA.InChiKey(predEPA.neighbors(i,:));
                res.EPA_DTXSID_neighbor(i,:)=train.CATMoS.model_EPA.DTXSID(predEPA.neighbors(i,:));
                %res.EPA_DSSTOXMPID_neighbor(i,:)=train.CATMoS.model_EPA.DSSTOXMPID(predEPA.neighbors(i,:));
                res.EPA_Exp_neighbor(i,:)=train.CATMoS.model_EPA.set.class_Exp(predEPA.neighbors(i,:));
                res.EPA_pred_neighbor(i,:)=train.CATMoS.model_EPA.set.class(predEPA.neighbors(i,:));
                
                res.GHS_CATMoS_ID_neighbor(i,:)=train.CATMoS.model_GHS.ChemID(predGHS.neighbors(i,:));
                res.GHS_CAS_neighbor(i,:)=train.CATMoS.model_GHS.CAS(predGHS.neighbors(i,:));
                res.GHS_InChiKey_neighbor(i,:)=train.CATMoS.model_GHS.InChiKey(predGHS.neighbors(i,:));
                res.GHS_DTXSID_neighbor(i,:)=train.CATMoS.model_GHS.DTXSID(predGHS.neighbors(i,:));
                %res.GHS_DSSTOXMPID_neighbor(i,:)=train.CATMoS.model_GHS.DSSTOXMPID(pred.neighbors(i,:));
                res.GHS_Exp_neighbor(i,:)=train.CATMoS.model_GHS.set.class_Exp(predGHS.neighbors(i,:));
                res.GHS_pred_neighbor(i,:)=train.CATMoS.model_GHS.set.class(predGHS.neighbors(i,:));
                
                res.LD50_CATMoS_ID_neighbor(i,:)=train.CATMoS.model_LD50.ChemID(predLD50.neighbors(i,:));
                res.LD50_CAS_neighbor(i,:)=train.CATMoS.model_LD50.CAS(predLD50.neighbors(i,:));
                res.LD50_InChiKey_neighbor(i,:)=train.CATMoS.model_LD50.InChiKey(predLD50.neighbors(i,:));
                res.LD50_DTXSID_neighbor(i,:)=train.CATMoS.model_LD50.DTXSID(predLD50.neighbors(i,:));
                %res.LD50_DSSTOXMPID_neighbor(i,:)=train.CATMoS.model_LD50.DSSTOXMPID(pred.neighbors(i,:));
                res.LD50_Exp_neighbor(i,:)=train.CATMoS.model_LD50.set.y_Exp(predLD50.neighbors(i,:));
                res.LD50_pred_neighbor(i,:)=train.CATMoS.model_LD50.set.y(predLD50.neighbors(i,:));
            end
            
            if strcmpi(ext,'.txt') && sep==1
                %res.Xtest=Xtest;
                fprintf(output(Locb(find(Locb))),'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output(Locb(find(Locb))),'VT experimental= %i, NT experimental= %i, EPA category experimental= %i, GHS category experimental= %i, LD50 experimental= %.3f\n', res.CATMoS_VT_exp(i),res.CATMoS_NT_exp(i),res.CATMoS_EPA_exp(i),res.CATMoS_GHS_exp(i),res.CATMoS_LD50_exp(i));
                end
                fprintf(output(Locb(find(Locb))),'VT predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_VT_pred(i),res.AD_VT(i),res.AD_index_VT(i),res.Conf_index_VT(i));
                fprintf(output(Locb(find(Locb))),'NT predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_NT_pred(i),res.AD_NT(i),res.AD_index_NT(i),res.Conf_index_NT(i));
                fprintf(output(Locb(find(Locb))),'EPA category predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_EPA_pred(i),res.AD_EPA(i),res.AD_index_EPA(i),res.Conf_index_EPA(i));
                fprintf(output(Locb(find(Locb))),'GHS category predicted= %i, AD: %i,AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_GHS_pred(i),res.AD_GHS(i),res.AD_index_GHS(i),res.Conf_index_GHS(i));
                fprintf(output(Locb(find(Locb))),'LD50 predicted= %.3f, AD: %i,AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_LD50_pred(i),res.AD_LD50(i),res.AD_index_LD50(i),res.Conf_index_LD50(i));
                
                if neighbors==1
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_VT.set.K, res.VT_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_NT.set.K, res.NT_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_EPA.set.K, res.EPA_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_GHS.set.K, res.GHS_CAS_neighbor{i,1:5});
                    fprintf(output(Locb(find(Locb))),'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_LD50.set.K, res.LD50_CAS_neighbor{i,1:5});
                end
                
                
            elseif strcmpi(ext,'.txt') && sep==0
                
                fprintf(output,'\t Molecule %s:\n', MoleculeNames{i});
                if exp
                    fprintf(output,'VT experimental= %i, NT experimental= %i, EPA category experimental= %i, GHS category experimental= %i, LD50 experimental= %.3f\n', res.CATMoS_VT_exp(i),res.CATMoS_NT_exp(i),res.CATMoS_EPA_exp(i),res.CATMoS_GHS_exp(i),res.CATMoS_LD50_exp(i));
                end
                fprintf(output,'VT predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_VT_pred(i),res.AD_VT(i),res.AD_index_VT(i),res.Conf_index_VT(i));
                fprintf(output,'NT predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_NT_pred(i),res.AD_NT(i),res.AD_index_NT(i),res.Conf_index_NT(i));
                fprintf(output,'EPA category predicted= %i, AD: %i, AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_EPA_pred(i),res.AD_EPA(i),res.AD_index_EPA(i),res.Conf_index_EPA(i));
                fprintf(output,'GHS category predicted= %i, AD: %i,AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_GHS_pred(i),res.AD_GHS(i),res.AD_index_GHS(i),res.Conf_index_GHS(i));
                fprintf(output,'LD50 predicted= %.3f, AD: %i,AD_index= %.2f, Conf_index= %.2f\n', res.CATMoS_LD50_pred(i),res.AD_LD50(i),res.AD_index_LD50(i),res.Conf_index_LD50(i));
                
                if neighbors==1
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_VT.set.K, res.VT_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_NT.set.K, res.NT_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_EPA.set.K, res.EPA_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_GHS.set.K, res.GHS_CAS_neighbor{i,1:5});
                    fprintf(output,'CAS of the %i nearest neighbors:%15s,%15s,%15s,%15s,%15s\n',train.CATMoS.model_LD50.set.K, res.LD50_CAS_neighbor{i,1:5});
                end
                
            end
        end
        
        
        if sep==1 && strcmpi(ext,'.csv')
            T=struct2table(res);
            if printtDesc==1
                %Xtest=[XtestVT; XtestNT; XtestEPA; XtestGHS; XtestLD50];
                Xtest=array2table(Xtest,'VariableNames',Desc);
                
                T=[T Xtest];
                res.Descriptors=Xtest;
            end
            writetable(T,FileOut{Locb(find(Locb))},'Delimiter',',');%,'QuoteStrings',true);
            fclose(output(Locb(find(Locb))));
            clear('T');
            
        elseif sep==0 && printtDesc==1 && strcmpi(ext,'.csv')
            
            %Xtest=[XtestVT; XtestNT; XtestEPA; XtestGHS; XtestLD50];
            Xtest(:,ismember(Desc,DescNames))=[];
            Desc(ismember(Desc,DescNames))=[];
            DescNames=[DescNames Desc];
            DescMat=[DescMat Xtest];
        end
        
        if sep==1
            resf.CATMoS=res;
            clear('res');
            
        end
        % Clean memory
        clear('XinCDK_CATMoS');
        clear('Xtest');
        clear('XtestVT');
        clear('XtestNT');
        clear('XtestEPA');
        clear('XtestGHS');
        clear('XtestLD50');
        clear('predVT');
        clear('predNT');
        clear('predEPA');
        clear('predGHS');
        clear('predLD50');
        clear('AD');
        %end clean memory
    end
    %--------------------------------------------------------------------------
    
    if sep==0 && strcmpi(ext,'.csv')
        res=struct2table(res);
        if printtDesc==1
            
            DescMat=array2table(DescMat,'VariableNames',DescNames);
            res=[res DescMat];
        end
        writetable(res,FileOut,'Delimiter',',');%,'QuoteStrings',true);
        fclose('all');
    end
    
    if sep==1
        res=resf;
        %res=0;
    end
    
    
    if verbose>0
        fprintf(1,'\n========== End Of Calculation ==========\n');
        fprintf('%i molecules predicted\n', length(MoleculeNames));
    end
    
    
end


