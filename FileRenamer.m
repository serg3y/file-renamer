function varargout = FileRenamer(varargin)

%Sergey Kharabash

%Compile:
% mbuild -setup                                                             %select compiler
% mcc -m    FileRenamer; delete mccExcludedFiles.log; delete readme.txt     %with dos shell 
% mcc -m -e FileRenamer; delete mccExcludedFiles.log; delete readme.txt     %no dos shell

%Issues:
% -If two files are to be renamed to same file name, only first will be reanemd.
% -Check for non allowed charecters before renaming
% -Use of ' in an operaor will break cfg file

%#ok<*DEFNU,*INUSD>
gui_Singleton = 1;
gui_State = struct(...
    'gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FileRenamer_OpeningFcn, ...
    'gui_OutputFcn',  @FileRenamer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function varargout = FileRenamer_OutputFcn(~,~,H)
varargout{1} = H.output;

% OPEN/CLOSE PROGRAM
function FileRenamer_OpeningFcn(h,~,H, varargin)
[F,f] = fileparts(mfilename('fullpath'));
diary(fullfile(F,[f '.log']))
diary on
fprintf('\nStarting: %s\n',datestr(now)) %#ok<TNOW1,DATST>
H.output = h;
try
    f = fopen('FileRenamer.ini','r');
    str(H.Folder,fgetl(f));
    str(H.Filter,fgetl(f));
    t = str2double(fgetl(f)); assert(any(t==0:1)); val(H.Recursive,t);
    t = str2double(fgetl(f)); assert(any(t==0:1)); val(H.Show_Extension,t);
    t = str2double(fgetl(f)); assert(any(t==0:1)); val(H.Show_Path,t);
    t = str2double(fgetl(f)); assert(any(t==1:4)); val(H.ChangeCase1,t);
    t = str2double(fgetl(f)); assert(any(t==1:4)); val(H.ChangeCase2,t);
    fclose(f);
catch e
    disp(e.message)
end
guidata(h, H);
t = findjobj(H.FindReplace);
set(H.FindReplace,'userd',t.getViewport.getView) %store java handle
UpdateOldList(H)

function FindReplace_CreateFcn(h,~,H)
try
    f = fopen('FileRenamer.cfg','r');
    t = eval(['{' fread(f,'*char')' '}']);
    t(:,1) = num2cell(cell2mat(t(:,1))>0); %convert 1|0 to true|false
    fclose(f);
catch e %use system default if file is currupt
    disp(e.message)
    t = repmat({false '' '' 'strrep'},[7 1]);
end
set(h,'data',t)

function FileRenamer_CloseRequestFcn(h,~,H)
[F,f] = fileparts(mfilename('fullpath'));
F = fullfile(F,f);
f = fopen([F '.ini'],'w');
fprintf(f,'%s\n',str(H.Folder));
fprintf(f,'%s\n',str(H.Filter));
fprintf(f,'%g\n',val(H.Recursive));
fprintf(f,'%g\n',val(H.Show_Extension));
fprintf(f,'%g\n',val(H.Show_Path));
fprintf(f,'%g\n',val(H.ChangeCase1));
fprintf(f,'%g\n',val(H.ChangeCase2));
fclose(f);
Save_Callback([],[],H,[F '.cfg'])
fprintf('Closed: %s\n\n',datestr(now)) %#ok<TNOW1,DATST>
diary off
delete(h);

%CALL BACKS
function Folder_Callback(~,~,H)
UpdateOldList(H)

function Open_Callback(~,~,H)
f = uigetdir(str(H.Folder),'Select Directory');
if ischar(f)
    str(H.Folder,f)
    UpdateOldList(H)
end

function Locate_Folder_Callback(~,~,H)
Navigate_Callback(H.Navigate,[],H)

function Up_one_Folder_Callback(~,~,H)
str(H.Folder,fileparts(str(H.Folder)))
UpdateOldList(H)

function Navigate_Callback(~,~,H)
system(['Explorer /e,' validfold(str(H.Folder)) ' &']);

function Filter_Callback(~,~,H)
UpdateOldList(H)

function Recursive_Callback(~,~,H)
UpdateOldList(H)

function Show_Extension_Callback(~,~,H)
UpdateOldList(H,val(H.List))

function Show_Path_Callback(~,~,H)
UpdateOldList(H,val(H.List))

function ChangeCase1_Callback(~,~,H)
UpdateNewList(H)

function ChangeCase2_Callback(~,~,H)
UpdateNewList(H)

function FindReplace_CellEditCallback(~,~,H)
UpdateNewList(H)

function List_Callback(h,~,H)
if numel(val(h)) == 1
    set(H.Manual,'en','on')
    NewList = get(H.List,'userd');
    [~,file,ext] = fileparts(NewList{val(h)});
    str(H.Manual,[file,ext])
else
    str(H.Manual,'','en','off')
end
UpdateNewList(H)

function Manual_Callback(h,~,H)
L = get(H.List,'userd'); %file list with full paths
file = L{val(H.List)};
RenameFile(file,fullfile(fileparts(file),str(h)))
UpdateOldList(H)

function Rename_Callback(~,~,H)
I = val(H.List); %selected source files
if ~isempty(I)
    L = get(H.List,'userd'); %source file names with full path and extension
    N = get(H.NewList,'userd'); %new file names
    II = find(~strcmp(L(I),N)); %skip if source and destination are same
%     t = unique([ALL_FILES;N(II)]);
%     if numel(unique(t))<numel(t) %check for confilcts
%         warning('dfg')
%         return
%     end
    for k = find(II)'
        RenameFile(L{I(k)},N{k})
    end
end
UpdateOldList(H)

% MENU CALLBACKS
function Locate_File_Callback(~,~,H)
f = get(H.List,'user');
if ~isempty(f)
    n = val(H.List);
    system(['explorer /select,' f{n(1)}]);
    %SHOpenFolderAndSelectItems can be used to select multiple files
end

function Save_Callback(~,~,H,f)
p = fullfile(fileparts(mfilename('fullpath')),'default');
if nargin<3
    [f,p] = uiputfile({'*.cfg' 'Config file (*.cfg)';'*.*' 'All files (*.*)'},'Save Operators',p);
    f = fullfile(p,f);
end
if ischar(f)
    d = get(H.FindReplace,'data')';
    d(2:3,:) = strrep(d(2:3,:),'''','''''');
    s = sprintf('%g\t''%s''\t''%s''\t''%s''\n',d{:});
    f = fopen(f,'wt');
    fwrite(f,s);
    fclose(f);
end

function Load_Callback(~,~,H)
p = fileparts(mfilename('fullpath'));
[f,p] = uigetfile({'*.cfg' 'Config file (*.cfg)';'*.*' 'All files (*.*)'},'Load Operators',p);
try
    f = fopen(fullfile(p,f),'r');
    t = fread(f,'*char')';
    fclose(f);
    data = eval(['{' t '}']);
    data(:,1) = cellfun(@(x)eq(x,1),data(:,1),'un',0); %double to logical
    set(H.FindReplace,'data',data)
catch e %file may be corrupt
    disp(e.identifier)
    disp(e.message)
end
UpdateNewList(H)

function Move_Up_Callback(~,~,H)
j = get(H.FindReplace,'userd');
l = j.getSelectedRows+1;
if ~isempty(l)
    j.clearSelection
    d = get(H.FindReplace,'data');
    I = 1:size(d,1);
    for k = 1:numel(l)
        if l(k)>1 && (k==1 || l(k-1)<l(k)-1)
            I([l(k) l(k)-1]) = I([l(k)-1 l(k)]);
            l(k) = l(k)-1;
        end
    end
    set(H.FindReplace,'data',d(I,:)); drawnow
    for k = 1:numel(l)
        j.changeSelection(l(k)-1,1,true,false);
        j.changeSelection(l(k)-1,2,true,true);
    end
end
UpdateNewList(H)

function Move_Down_Callback(~,~,H)
j = get(H.FindReplace,'userd');
l = j.getSelectedRows+1;
if ~isempty(l)
    j.clearSelection
    d = get(H.FindReplace,'data');
    I = 1:size(d,1);
    for k = numel(l):-1:1
        if l(k)<numel(I) && (k==numel(l) || l(k+1)>l(k)+1)
            I([l(k) l(k)+1]) = I([l(k)+1 l(k)]);
            l(k) = l(k)+1;
        end
    end
    set(H.FindReplace,'data',d(I,:)); drawnow
    for k = 1:numel(l)
        j.changeSelection(l(k)-1,1,true,false);
        j.changeSelection(l(k)-1,2,true,true);
    end
end
UpdateNewList(H)

function Insert_Callback(~,~,H)
d = get(H.FindReplace,'data');
j = get(H.FindReplace,'userd'); %selected cells, n*m, can be empty
i = j.getSelectedRows+1;
if isempty(i)
    n = 1; %add 1
    i = size(d,1); %at the end
else
    n = numel(i(:,1)); %add n
    i = min(i(:,1))-1; %before selection
end
d = [d(1:i,:);repmat({false '' '' 'strrep'},n,1);d(i+1:end,:)];
set(H.FindReplace,'data',d); drawnow
for i = i:i+n-1
    j.changeSelection(i,1,true,false);
    j.changeSelection(i,2,true,true);
end

function Delete_Callback(~,~,H)
j = get(H.FindReplace,'userd');
i = j.getSelectedRows+1;
if ~isempty(i)
    d = get(H.FindReplace,'data');
    d(unique(i(:,1)),:) = [];
    set(H.FindReplace,'data',d);
end
UpdateNewList(H)

function Help_Operators_Callback(h,e,H) %under construction
dialog('name','Example Operators');

% MAIN FUNCTIONS
function UpdateOldList(H,I)
fold = str(H.Folder);
filt = str(H.Filter);
if ~isfolder(fold)
    set(H.Folder,'back',[1 .8 .8]) %red, invalid path
    val(H.List,[],'userd',{})
else
    set(H.Folder,'back',[1 1 1]) %white, good path
    paths = FileSearchR(fold,filt,val(H.Recursive)); %list of file paths before renaming
    [F,files,e] = cellfun(@fileparts,paths,'un',0);
    if val(H.Show_Extension) && ~isempty(files)
        files = strcat(files,e);
    end
    if val(H.Show_Path)
        F = strrep(F,fold,''); %show relative path, HACK!
        files = cellfun(@fullfile,F,files,'un',0);
    end
    str(H.List,files,'userd',paths) %will need the full path when renaming
    if nargin<2 || isempty(I)
        I = 1:numel(files); %select all
    end
    val(H.List,I) %select entries
end
str(H.Manual,'','enable','off') %click a file to turn on
set(H.FileRenamer,'Name',[fold ' - FileRenamer']) %program title
UpdateNewList(H)
function UpdateNewList(H)
I = val(H.List);
if isempty(I)
    str(H.NewList,{},'userd',{}) %clear display & data
    val(H.List,[]) %clear selection
else
    d = get(H.List,'userd');
    [F,files,e] = cellfun(@fileparts,d(I),'un',0);
    if val(H.Show_Extension)
        files = strcat(files,e);
    end
    files = ChangeCase(files,val(H.ChangeCase1));
    rules = get(H.FindReplace,'data'); %find replace rules
    for k = find([rules{:,1}])
        switch rules{k,4}
            case 'strrep',   files = strrep(files,rules{k,2},rules{k,3});
            case 'regexp',   files = regexp(files,rules{k,2},'match','once'); %fix needed
            case 'regexprep',files = regexprep(files,rules{k,2},rules{k,3}); %fix needed
            case 'prepend',  files = strcat(rules{k,3},files);
            case 'append',   files = strcat(files,rules{k,3});
        end
    end
    files = ChangeCase(files,val(H.ChangeCase2));
    str(H.NewList,files) %display
    if val(H.Show_Extension)
        paths = cellfun(@fullfile,F,files,'un',0);
    else
        paths = cellfun(@fullfile,F,strcat(files,e),'un',0);
    end
    val(H.NewList,[],'userd',paths) %selection & data
end
str(H.Info,sprintf('%g files',numel(val(H.NewList)))) %number of files selected for renaming

% FUNCTION HELPERS
function I = shift(l,d)
%Generates shifted indesies based on a selection, eg shift([1 0 0 1],-1)
if nargin<2 || isempty(d), d = 1; end
if d>0
    I = 1:numel(l);
    for k = 1:numel(I)-1
        if ~l(k) && l(k+1)
            I([k k+1]) = I([k+1 k]);
            l([k k+1]) = [1 0];
        end
    end
else
    I = 1:numel(l);
    for k = numel(I):-1:2
        if ~l(k) && l(k-1)
            I([k k-1]) = I([k-1 k]);
            l([k k-1]) = [1 0];
        end
    end
end
function fold = validfold(fold)
while ~isfolder(fold)
    fold = fileparts(fold);
    if isempty(fold)
        fold = 'C:\'; break
    end
end

% Folder = 'C:\'; %can be a relative path
% Names_Only = arrayfun(@char,java.io.File(Folder).list,'un',0) %cellstr
% Full_Paths = arrayfun(@char,java.io.File(Folder).listFiles,'un',0) %cellstr
% 
% %list files only in a folder (skip folders)
% Folder = 'C:\';
% jPaths = java.io.File(Folder).listFiles; %java.io.File objects
% jNames = java.io.File(Folder).list; %java.lang.String objects
% isFolder = arrayfun(@isDirectory,jPaths); %boolean
% File_Names_Only = arrayfun(@char,jNames(~isFolder),'un',0) %cellstr


function F = FileSearchR(fldr,wild,recursive)
F = {};
D = dir(fldr);
L = strcmp({D(:).name},'.') | strcmp({D(:).name},'..'); %ignire . & ..
if recursive
    t = sort_nat({D([D.isdir] & ~L).name});
    for k = 1:numel(t) %traverse sub folders
        %button = questdlg('qstring',fldr,'Stop','Stop'); %under construction, would like to be able to abort if the list is too large
        F = [F; FileSearchR(fullfile(fldr,t{k}),wild,1)]; %#ok<AGROW>
    end
end
wild = strrep(strrep(strrep(wild,',','|'),'.','\.'),'*','.*'); %wildcard as regexp, probably could have used regexptranslate('wildcard', '*.bmp') somehow
t = {D(~[D.isdir] &~ L).name}'; %all files in folder
L = strcmp(t,regexp(t,wild,'match','once')); %matches
if ~any(L)
    return
end
if strcmp(fldr(end),filesep)
    fldr = fldr(1:end-1); %avoid double slashes 'fodler\\file'
end
F = [F; strcat(fldr,filesep,sort_nat(t(L)))]; %catinate matching files
function str = ChangeCase(str,type)
switch lower(type(1))
    case {'l' 2} %lower
        str = lower(str);
    case {'t' 3} %Title
        if iscell(str)
            str = cellfun(@(x)ChangeCase(x,3),str,'un',0); return %hack
        end
        str = lower(str); %lower the case
        a = isstrprop(str,'alpha'); %letter
        s = find(diff([0 a]) > 0); %word start
        if numel(s)>0
            e = find(diff([a 0]) < 0); %word end
            w = mat2cell(str(a),1,e-s+1); %words
            l = {'at' 'the' 'and' 'on' 'to' 'in' 'with' 'but' 'or' 'is'}; %exceptions
            t = cellfun(@(x)strcmp(w(2:end),x),l,'un',0); %compare with short words
            i = [0 any(cat(1,t{:}),1)]; %confirmed short words
            str(s(~i)) = upper(str(s(~i)));
        end
    case {'u' 4} %UPPER
        str = upper(str);
end
function o = val(h,varargin)
if nargout
    o = get(h,'val');
end
if nargin>1
    set(h,'val',varargin{:});
end
function o = str(h,varargin)
if nargout
    o = get(h,'str');
end
if nargin>1
    set(h,'str',varargin{:});
end
function RenameFile(f1,f2)
if ~strcmp(f1,f2)
    [F,t1,e1] = fileparts(f1);
    [~,t2,e2] = fileparts(f2);
    if java.io.File(f1).renameTo(java.io.File(f2))
        fprintf('"%s" > "%s" : "%s" \n',[t1 e1],[t2 e2],F)
    else
        fprintf('FAILED: "%s" > "%s" : "%s" \n',[t1 e1],[t2 e2],F)
    end
end
function list = jls(fold,type,mode)
%java based ls equivilent

%defaults
if nargin<1 || isempty(fold), fold = cd; end
if nargin<2 || isempty(type), type = ''; end
if nargin<3 || isempty(mode), mode = 'full'; end

%get list
jfold = java.io.File(fold);
jlist = jfold.list;

%filter
if ~isempty(type)
    B = arrayfun(@(x)x.endsWith(type),jlist); %boolean
    jlist = jlist(B);
end

%convert to cellstr
list = arrayfun(@char,jlist,'un',0);

%full paths
switch mode
    case 'short'
        %nothing t do
    case 'full'
        list = strcat(char(jfold),filesep,list);
end