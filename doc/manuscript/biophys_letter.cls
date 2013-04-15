\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{biophys_letter}
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\fi
\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}



\DeclareOption{cmfonts}{%
    \g@addto@macro{\@greeksetup}{%
%        \DeclareMathSymbol\varGamma    {\mathord}{letters}{0}
%        \DeclareMathSymbol\varDelta    {\mathord}{letters}{1}
%        \DeclareMathSymbol\varTheta    {\mathord}{letters}{2}
%        \DeclareMathSymbol\varLambda   {\mathord}{letters}{3}
%        \DeclareMathSymbol\varXi       {\mathord}{letters}{4}
%        \DeclareMathSymbol\varPi       {\mathord}{letters}{5}
%        \DeclareMathSymbol\varSigma    {\mathord}{letters}{6}
%        \DeclareMathSymbol\varUpsilon  {\mathord}{letters}{7}
%        \DeclareMathSymbol\varPhi      {\mathord}{letters}{8}
%        \DeclareMathSymbol\varPsi      {\mathord}{letters}{9}
%        \DeclareMathSymbol\varOmega    {\mathord}{letters}{10}
    }
}

\DeclareOption{slantedgreek}{\def\@greeksetup{\AtBeginDocument{%
% \renewcommand{\Gamma}{{\mit\varGamma}}
% \renewcommand{\Delta}{{\mit\varDelta}}
% \renewcommand{\Theta}{{\mit\varTheta}}
% \renewcommand{\Lambda}{{\mit\varLambda}}
% \renewcommand{\Xi}{{\mit\varXi}}
% \renewcommand{\Pi}{{\mit\varPi}}
% \renewcommand{\Sigma}{{\mit\varSigma}}
% \renewcommand{\Upsilon}{{\mit\varUpsilon}}
% \renewcommand{\Phi}{{\mit\varPhi}}
% \renewcommand{\Psi}{{\mit\varPsi}}
% \renewcommand{\Omega}{{\mit\varOmega}}
}}}
\DeclareOption{uprightgreek}{\def\@greeksetup{}}

\ExecuteOptions{letterpaper,10pt,twoside,twocolumn,final,slantedgreek,cmfonts}
\ProcessOptions
%\input{size1\@ptsize.clo}

%\JOURNAL@fonts
\RequirePackage{graphicx}
\RequirePackage{latexsym,amsmath,lastpage,amsthm,fancybox}
\usepackage{times}

% \DeclareMathAlphabet\mathbi
%   \encodingdefault\rmdefault\bfdefault\itdefault
% \DeclareMathAlphabet\mathsc
%   \encodingdefault{cmr}\mddefault\scdefault

\renewcommand\normalsize{%
  \@setfontsize\normalsize{10}{12}%
  \abovedisplayskip 10\p@ \@plus 2\p@ \@minus5\p@
  \abovedisplayshortskip \z@ \@plus 3\p@
  \belowdisplayshortskip 6\p@ \@plus 3\p@ \@minus3\p@
 \parindent 10pt \belowdisplayskip \abovedisplayskip
  \let\@listi\@listI\baselineskip12pt
  }

\normalsize
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\ack{\vspace*{2pt}
   \@setfontsize\small\@viiipt{10}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\ssmall{%
   \@setfontsize\ssmall\@ixpt{10}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{9.5}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}


% \newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
% \newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
% \newcommand\large{\@setfontsize\large\@xipt{13.5}}
% \newcommand\Large{\@setfontsize\Large\@xivpt{18}}
% \newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
% \newcommand\huge{\@setfontsize\huge\@xxpt{25}}
% \newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0pt\relax}
\setlength\parindent{12\p@}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@twosidetrue
\@mparswitchtrue
% \setlength\headheight{12\p@}
% \setlength\headsep   {0.25in}
% \setlength\topskip   {10\p@}
% \setlength\footskip{30\p@}
% \setlength\maxdepth{0.5\topskip}
% \setlength\textwidth{35pc}
%   \@settopoint\textwidth
% \setlength\textheight{47\baselineskip}
% \addtolength\textheight{\topskip}
%   \@settopoint\textheight
% \setlength\marginparsep{11\p@}
% \setlength\marginparpush{5\p@}
% \setlength\oddsidemargin   {44\p@}
% \setlength\evensidemargin  {82\p@}
% \setlength\marginparwidth {107\p@}
%     \setlength\@tempdima        {\paperwidth}
%     \addtolength\@tempdima      {-\textwidth}
%     \setlength\oddsidemargin    {.5\@tempdima}
%     \addtolength\oddsidemargin  {-1in}
%     \setlength\marginparwidth   {.5\@tempdima}
%     \addtolength\marginparwidth {-\marginparsep}
%     \addtolength\marginparwidth {-0.4in}
%     \addtolength\marginparwidth {-.4in}
%   \ifdim \marginparwidth >2in
%      \setlength\marginparwidth{2in}
%   \fi
%   \@settopoint\oddsidemargin
%   \@settopoint\marginparwidth
%   \setlength\evensidemargin  {\paperwidth}
%   \addtolength\evensidemargin{-2in}
%   \addtolength\evensidemargin{-\textwidth}
%   \addtolength\evensidemargin{-\oddsidemargin}
%   \@settopoint\evensidemargin
%   \setlength\topmargin{\paperheight}
%   \addtolength\topmargin{-2in}
%   \addtolength\topmargin{-\headheight}
%   \addtolength\topmargin{-\headsep}
%   \addtolength\topmargin{-\textheight}
%   \addtolength\topmargin{-\footskip}     % this might be wrong!
%   \addtolength\topmargin{-0.5\topmargin}
%   \@settopoint\topmargin


%%%%%%%%%%%%%%%%%%%%trim size setting

\columnsep 18pt
\textwidth41.6pc
\textheight54.85pc
\oddsidemargin-1.5pc
\evensidemargin-1.7pc
\topmargin4.1mm
\headsep20pt
\parindent15pt
\footskip20pt 

\@greeksetup


\def\gridframe#1{
\let\grid=#1
\ifx Y\grid
\def\gridpl{
\setlength{\fboxsep}{0pt}\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(-0.14,1.81){\framebox(41.52,54.35)}%0.1,2%-1,3
\end{picture}
}
\def\gridod{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(-0.01,-56){\framebox(41.57,54.55)}
\end{picture}}
\def\gridev{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(-0.01,-56){\framebox(41.6,54.55)}
\end{picture}}
\else
\def\gridpl{}
\def\gridod{}
\def\gridev{}
\fi
}


\def\ob{{\rm [}}
\def\cb{{\rm ]}}

%%%%%%%%%%%%%%%Plain Style

\def\jno#1{\gdef\@kxi{#1}}


% \def\copyline#1{
% \let\copyl=#1
% \ifx P\copyl
% \def\copyl{\begin{tabular}{@{}lr@{}}
% \noindent\sffamily\fontsize{6.8}{9}\selectfont Biophysical Journal: Biophysical
% Letters
% \end{tabular}\hfill
% }
% \else
% \def\copyl{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr}\\[-14pt]
% \hspace*{-3pt}\sffamily\fontsize{6.8}{9}\selectfont 
% Biophysical Journal: Biophysical Letters & \noindent\sffamily\fontsize{6.8}{9}\selectfont %L0{\thepage}
% \end{tabular*}\hfill
% }
% \fi
% }



%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\jno#1{\gdef\@kxi{#1}}


\def\ps@plain{
\def\@oddfoot{\gridpl{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[-6pt]
\sffamily\fontsize{8}{10}\selectfont 
Biophysical Journal: Biophysical Letters & \noindent\sffamily\fontsize{8}{10}\selectfont L0{\thepage}
\end{tabular*}\hfill
}}
      \def\@evnhead{\CROPMARKSA}
           \def\@oddhead{\CROPMARKSA \begin{tabular}{@{\hspace*{-2.3pt}}lr@{}}\\[-14pt]
           \fontsize{8}{10}\sffamily\selectfont{
           Biophysical Journal: Biophysical Letters}
\end{tabular}}
       \let\@mkboth\markboth
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
      




\def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
     \def\@evenfoot{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[-6pt]
     {\fontsize{8}{10}\sffamily\selectfont\leftmark} & {\fontsize{8}{10}\sffamily\selectfont{L0\thepage}}
     \end{tabular*}}
     \def\@oddfoot{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[-6pt]
     {\fontsize{8}{10}\sffamily\selectfont\leftmark} & {\fontsize{8}{10}\sffamily\selectfont{L0\thepage}}
     \end{tabular*}}
      \def\@evenhead{\CROPMARKS\gridev{\fontsize{8}{10}\sffamily\selectfont\leftmark}\hfill}
      \def\@oddhead{\CROPMARKS\gridod{\fontsize{8}{10}\sffamily\selectfont\leftmark}\hfill}
      \let\@mkboth\markboth
        }

\setlength\footnotesep{6.65\p@}
\setlength{\skip\footins}{0.6pc}%{1pc}%krishna
\setcounter{topnumber}{8}
\renewcommand\topfraction{1}
\setcounter{bottomnumber}{3}
\renewcommand\bottomfraction{0.8}
\setcounter{totalnumber}{8}
\renewcommand\textfraction{0}
\renewcommand\floatpagefraction{0.8}
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\@fptop{0\p@ \@plus 1fil}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\let\@accepted\@empty
\def\accepted#1{\gdef\@accepted{#1}}
\let\@received\@empty
\def\received#1{\gdef\@received{#1}}
\let\@revised\@empty
\def\revised#1{\gdef\@revised{#1}}
\def\volume#1{\setcounter{@volumenumber}{#1}}
  \newcounter{@volumenumber}
  \setcounter{@volumenumber}{1}
\def\issue#1{\setcounter{@issuenumber}{#1}}
  \newcounter{@issuenumber}
  \setcounter{@issuenumber}{1}
\def\journalname#1{\gdef\@journalname{#1}}
\def\journalnumber#1{\gdef\@journalnumber{#1}}
\def\copyrightname#1{\gdef\@copyrightname{#1}}
\def\price#1{\gdef\@price{#1}}
\volume{1}
\issue{1}
\journalname{Biostatistics}
\journalnumber{0006--3444}
\copyrightname{Oxford University Press}
\price{\$124.00}
\def\label@firstpage{%
   \immediate\write\@auxout{\string
   \newlabel{firstpage}{{}{\thepage}}}}
\def\label@lastpage{%
  \advance\c@page\m@ne
   \immediate\write\@auxout{\string
   \newlabel{lastpage}{{}{\thepage}}}\advance\c@page\@ne}
  \newlabel{lastpage}{{}{0}}
\def\ps@biostatistics{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty%
      \def\@evenhead{{\thepage\hfil\itshape
      {\@shortauthors\@andorcomma\relax\@andorcomma\relax}\relax}\hfil}%
      \def\@oddhead{\hfil{\itshape\@shorttitle}\hfil\thepage}
     \relax
     \relax}%
\def\testpages{\expandafter\@secondoftwo\r@firstpage}
\def\pagestest{\expandafter\@secondoftwo\r@lastpage}
\def\ps@title{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty%
      \def\@oddhead{\CROPMARKSA%
      \smash{\lower12pt\rlap{\small\itshape Printed in Great Britain}}%
       \vbox{\hsize\textwidth\@parboxrestore\small
         {\itshape\@journalname} ({\the\year}),
          \textbf{\arabic{@volumenumber}},
          \arabic{@issuenumber}, \textit{pp}.~\pageref{firstpage}%
         \ifnum\testpages=\pagestest\else
        --\pageref{lastpage}\fi\hfil}
     \relax
     \relax}%
     \let\@evenhead\@oddhead}
\def\ps@empty{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \let\@evenhead@empty\let\@oddhead\@empty
  }
\newif\if@@address
\@@addressfalse
\let\@authors\@empty
\def\address#1{%
  \ifx\@authors\@empty
    \def\@authors{\sffamily{#1}}
  \else
    \g@addto@macro\@authors{\if@@address\\*\else\\*[4\p@]\fi
    {\begin{tabular}{@{}l@{}}
    \parbox{175mm}{\vbox{\fontsize{8}{9}\sffamily#1}}
    \end{tabular}}}
  \fi}
\def\affiliation#1{#1}
\def\email#1{%
  \ifx\@authors\@empty
    \def\@authors{{#1}{\@@addresstrue}}
  \else
    \g@addto@macro\@authors{\\* {#1}\@@addresstrue}
  \fi}
\def\title{\@dblarg\JOURNAL@title}
\def\JOURNAL@title[#1]#2{%
  \def\@title{{#2}}%
  \def\@shorttitle{#1}}
\newsavebox{\@and}
\savebox{\@and}{{AND}}
\def\and{\usebox{\@and}}
\def\author{\@dblarg\JOURNAL@author}
\def\JOURNAL@author[#1]#2{%
  \ifx\@authors\@empty
    \def\@authors{{\begin{tabular}{@{}l@{}}
    \parbox{175mm}{\vbox{\fontsize{10}{12}\sffamily#1}}
    \end{tabular}}\\[-13.5pt]}%
  \else
    \g@addto@macro\@authors{\if@@address\\*[8pt]\else\\*\fi\@@addressfalse}%
    \g@addto@macro\@authors{{\MakeUppercase{#2}}}%
  \fi\biostatistics@addshortauthor{#1}}
\def\biostatistics@addshortauthor#1{%
  \g@addto@macro\@shortauthors{%
    \@andorcomma,{\ignorespaces#1\unskip}}}
\def\authorsetal{%
  \let\biostatistics@addshortauthor\@gobble
  \g@addto@macro\@shortauthors{%
    \@andorcomma,\etalname}}
\def\etalname{\textit{et~al}.}
\let\@authors\@empty
\def\@shortauthors{\relax\@gobbletwo}
\def\@andorcomma#1#2\@andorcomma#3{%
  \ifx#3\relax
    \ifx#1\relax%
      \let\@tempa\relax
    \else
      \def\@tempa{#2}\def\@tempb{\etalname}%
      \ifx\@tempa\@tempb
        \def\@tempa{\space#2\@gobbletwo}%
      \else
        \def\@tempa{\space\andname\space #2\@gobbletwo}%
      \fi
    \fi
  \else
    \def\@tempa{#1 #2\@andorcomma#3}%
  \fi
  \@tempa}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%################%%%%%%%%%%%

%%%%%%%%%%#############%%%%%%%%%%%%%%
\newcommand\maketitle{%
 \thispagestyle{plain}\label@firstpage
  \renewcommand{\thefootnote}{\fnsymbol{footnote}}
    \setcounter{footnote}{0}
  \let\footnoterule\relax
  \let \footnote \thanks
  \null
  \vbox{\onecolumn
    {\vspace*{-15pt}\titlefont\@title\par}%
    \vspace*{2pc}%
    \noindent%
    {\begin{tabular}[t]{@{}l}
        \fontsize{10}{12}\@authors
      \end{tabular}\par}%
      %\end{center}
\@thanks
  \null
   \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@authors\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
  \vspace*{0pt}
  }\enlargethispage{-12pt}
}


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%####################
\setcounter{secnumdepth}{3}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\renewcommand\thesection      {\@arabic\c@section}
\renewcommand\thesubsection   {{\rm \thesection.\@arabic\c@subsection}}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection\@arabic\c@paragraph}
%\def\@seccntformat#1{\csname the#1\endcsname\enspace }

\def\titlefont{\fontsize{14}{15}\bfseries\sffamily\noindent}







%\newcount\@secpenalty
\@secpenalty = -300
 
 
\newif\if@noskipsec \@noskipsectrue


\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
 
\def\@newstartsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@newsect{#1}{#2}{#3}{#4}{#5}{#6}}}}
 
 

\def\@newsect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname.\hskip4pt}\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M  #8\par}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname}\fi
                       #7}}\fi
     \@xsect{#5}}
 
\def\@xsect#1{\@tempskipa #1\relax
      \ifdim \@tempskipa>\z@
       \par \nobreak
       \vskip \@tempskipa
       \@afterheading
    \else \global\@nobreakfalse \global\@noskipsectrue
       \everypar{\if@noskipsec \global\@noskipsecfalse
                   \clubpenalty\@M \hskip -\parindent
                   \begingroup \@svsechd \endgroup \unskip
                   \hskip -#1\relax  % relax added 14 Jan 91
                  \else \clubpenalty \@clubpenalty
                    \everypar{}\fi}\fi\ignorespaces}
 
% \def\@ssect#1#2#3#4#5{\@tempskipa #3\relax
%    \ifdim \@tempskipa>\z@
%      \begingroup #4\@hangfrom{\hskip #1}{\interlinepenalty \@M #5\par}\endgroup
%    \else \def\@svsechd{#4\hskip #1\relax #5}\fi
%     \@xsect{#3}}
%  
\newif\if@afterindent \@afterindenttrue
 
\def\@afterheading{\global\@nobreaktrue
      \everypar{\if@nobreak
                   \global\@nobreakfalse
                   \clubpenalty \@M
                   \if@afterindent \else {\setbox\z@\lastbox}\fi
                 \else \clubpenalty \@clubpenalty
                    \everypar{}\fi}}
 
 
% \@hangfrom{TEXT} : Puts TEXT in a box, and makes a hanging indentation
%    of the following material up to the first \par.  Should be used
%    in vertical mode.
%
\def\@hangfrom#1{\setbox\@tempboxa\hbox{#1}%
      \hangindent \wd\@tempboxa\noindent\box\@tempboxa}
     
     
     
     
     



\newcommand\section{\@newstartsection {section}{1}{\z@}%
                                   {1.75pc\@plus0.125\p@\@minus0.125\p@}%
                                   {5.5pt\relax}%
                                   {\normalfont\normalsize\bfseries\sffamily}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {1.75pc\@plus0.125\p@\@minus0.125\p@}%
                                     {6pt\relax}%
                                     {\normalfont\itshape\centering\noindent}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {1.25pc\@plus0.0125\p@\@minus0.0125\p@}%
                                     {0.25pc \@plus0.0125\p@}%
                                     {\normalfont\normalsize\itshape\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                     {1pc\@plus0.0125\p@\@minus0.0125\p@}%
                                     {-8.5\p@}%
                                     {\normalfont\normalsize\itshape}}
  \setlength\leftmargini  {2.5em}
\leftmargin  \leftmargini
\setlength\leftmarginii  {24\p@}
\setlength\leftmarginiii {18\p@}
\setlength\leftmarginiv  {16\p@}
\setlength\leftmarginv   {12\p@}
\setlength\leftmarginvi  {12\p@}
\setlength  \labelsep  {0.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\setlength\partopsep{2\p@ \@plus 1\p@ \@minus 1\p@}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\def\@listi{\leftmargin.5\leftmargini
            \parsep 4\p@ \@plus2\p@ \@minus\p@
            \topsep 8\p@ \@plus2\p@ \@minus4\p@
            \itemsep-3\p@ \@plus2\p@ \@minus\p@}
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    4\p@ \@plus2\p@ \@minus\p@
              \parsep    2\p@ \@plus\p@  \@minus\p@
              \itemsep   \parsep}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    2\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep   \topsep}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}

\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

\newenvironment{abstract}[3]{\vbox{\fontsize{9}{11}\selectfont\vspace*{16pt}\noindent{{\fontsize{9}{11}\sffamily\selectfont
ABSTRACT}\quad}{\fontsize{9}{10}\sffamily#1}\\[15pt]
     {\centerline{\fontsize{10}{11}\sffamily\selectfont{#2}}}\\[9pt]
     {\centerline{\fontsize{9}{11}\sffamily\selectfont{#3}}}}}{
}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}

\newbox\boxH
\newbox\boxHrunin


\newenvironment{NL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\arabic{enumii}.}{%
 \usecounter{enumii}
  \topsep=6pt %plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=24pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%



%%%%%%%%%%%%%%%%%%%%bracketNL
%
%

\newenvironment{BNL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\arabic{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=6pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=24pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%



%%%%%%%%%%%%%%%%%AlphabetNL
%
%


\newenvironment{ANL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\alph{enumii}.}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=12pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%




%%%%%%%%%%%%%%roman NL
%
%

\newenvironment{RNL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\roman{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=24pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%


%%%%%%%%%%%DRNL
%
\newenvironment{DRNL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{(\roman{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=28pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%



%%%%%%%%%%Double alpha NL
%
%
\newenvironment{DANL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{(\Alph{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=32pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%



\newenvironment{BL}{
  \begin{list}{$\bullet$}{\fontsize{10}{12}\selectfont%
 \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=6pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\labelsep
 \leftmargin=12pt
 \rightmargin=0pt
 \listparindent=0pt}
}
{\end{list}}%




\newenvironment{UNL}{
  \begin{list}{}{
 \topsep=0pt %6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 %\labelsep=5pt
 %\labelwidth\labelsep
 \leftmargin=0pt
 \rightmargin=0pt
 \listparindent=0pt}
}
{\end{list}}%




% \def\newdefinition#1{%
%   \@ifnextchar[{\@odef{#1}}{\@ndef{#1}}}
% \def\@ndef#1#2{%
%   \@ifnextchar[{\@xndef{#1}{#2}}{\@yndef{#1}{#2}}}
% \def\@xndef#1#2[#3]{%
%   \expandafter\@ifdefinable\csname #1\endcsname
%     {\@definecounter{#1}\@newctr{#1}[#3]%
%      \expandafter\xdef\csname the#1\endcsname{%
%        \expandafter\noexpand\csname the#3\endcsname \@defcountersep
%           \@defcounter{#1}}%
%      \global\@namedef{#1}{\@def{#1}{#2}}%
%      \global\@namedef{end#1}{\@enddefinition}}}
% \def\@yndef#1#2{%
%   \expandafter\@ifdefinable\csname #1\endcsname
%     {\@definecounter{#1}%
%      \expandafter\xdef\csname the#1\endcsname{\@defcounter{#1}}%
%      \global\@namedef{#1}{\@def{#1}{#2}}%
%      \global\@namedef{end#1}{\@enddefinition}}}
% \def\@odef#1[#2]#3{%
%   \@ifundefined{c@#2}{\@nocounterr{#2}}%
%     {\expandafter\@ifdefinable\csname #1\endcsname
%     {\global\@namedef{the#1}{\@nameuse{the#2}}%
%   \global\@namedef{#1}{\@def{#2}{#3}}%
%   \global\@namedef{end#1}{\@enddefinition}}}}
% \def\@def#1#2{%
%   \refstepcounter{#1}%
%   \@ifnextchar[{\@ydef{#1}{#2}}{\@xdef{#1}{#2}}}
% \def\@xdef#1#2{%
%   \@begindefinition{#2}{\csname the#1\endcsname}\ignorespaces}
% \def\@ydef#1#2[#3]{%
%   \@opargbegindefinition{#2}{\csname the#1\endcsname}{#3}\ignorespaces}
% \def\@defcounter#1{\noexpand\arabic{#1}}
% \def\@defcountersep{.}
% \def\@begindefinition#1#2{\trivlist
%    \item[\hskip \labelsep{\indent\scshape #1\ #2.}]\rmfamily}
% \def\@opargbegindefinition#1#2#3{\trivlist
%       \item[\hskip \labelsep{\indent\scshape #1\ #2\ (#3).}]\rmfamily}
% \def\@enddefinition{\endtrivlist}
% \def\@begintheorem#1#2{\trivlist
%    \item[\hskip \labelsep{\indent\scshape #1\ #2.}]\itshape}
% \def\@opargbegintheorem#1#2#3{\trivlist
%       \item[\hskip \labelsep{\indent\scshape #1\ #2\ (#3).}]\itshape}
% \def\newproof#1{%
%   \@ifnextchar[{\@oprf{#1}}{\@nprf{#1}}}
% \def\@nprf#1#2{%
%   \@ifnextchar[{\@xnprf{#1}{#2}}{\@ynprf{#1}{#2}}}
% \def\@xnprf#1#2[#3]{%
%   \expandafter\@ifdefinable\csname #1\endcsname
%     {\@definecounter{#1}\@newctr{#1}[#3]%
%      \expandafter\xdef\csname the#1\endcsname{%
%        \expandafter\noexpand\csname the#3\endcsname \@prfcountersep
%           \@prfcounter{#1}}%
%      \global\@namedef{#1}{\@prf{#1}{#2}}%
%      \global\@namedef{end#1}{\@endproof}}}
% \def\@ynprf#1#2{%
%   \expandafter\@ifdefinable\csname #1\endcsname
%     {\@definecounter{#1}%
%      \expandafter\xdef\csname the#1\endcsname{\@prfcounter{#1}}%
%      \global\@namedef{#1}{\@prf{#1}{#2}}%
%      \global\@namedef{end#1}{\@endproof}}}
% \def\@oprf#1[#2]#3{%
%   \@ifundefined{c@#2}{\@nocounterr{#2}}%
%     {\expandafter\@ifdefinable\csname #1\endcsname
%     {\global\@namedef{the#1}{\@nameuse{the#2}}%
%   \global\@namedef{#1}{\@prf{#2}{#3}}%
%   \global\@namedef{end#1}{\@endproof}}}}
% \def\@prf#1#2{%
%   \refstepcounter{#1}%
%   \@ifnextchar[{\@yprf{#1}{#2}}{\@xprf{#1}{#2}}}
% \def\@xprf#1#2{%
%   \@beginproof{#2}{\csname the#1\endcsname}\ignorespaces}
% \def\@yprf#1#2[#3]{%
%   \@opargbeginproof{#2}{\csname the#1\endcsname}{#3}\ignorespaces}
% \def\@prfcounter#1{\noexpand\arabic{#1}}
% \def\@prfcountersep{.}
% \def\@beginproof#1#2{\trivlist
%    \item[\hskip \labelsep{\indent\scshape #1.}]\rmfamily}
% \def\@opargbeginproof#1#2#3{\trivlist
%       \item[\hskip \labelsep{\indent\scshape #1\ (#3).}]\rmfamily}
% \def\@endproof{\endtrivlist}
% \newcommand*{\qed}{\hbox{}\hfill$\Box$}
% \newproof{proof}{Proof}
\newcommand\appendix{
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \renewcommand\thesection{\@Alph\c@section}%
  \setcounter{equation}{0}%
  \@addtoreset{equation}{section}%
\def\@Alph##1{%
  \ifcase##1 A\or A\or B\or C\or D\or E\or F\or G\or H\or I\or J\or
   K\or L\or M\or N\or O\or P\or Q\or R\or S\or T\or U\or V\or W\or X\or
   Y\or Z\else\@ctrerr\fi}%
  \renewcommand\theequation{\@Alph\c@section.\@arabic\c@equation}%
\section*{APPENDIX}}

\newcommand\appendixone{
%
  %\renewcommand\thesection{\@Alph\c@section.}%
  \setcounter{equation}{0}%
  \@addtoreset{equation}{section}%
\def\@Alph##1{%
  \ifcase##1 A\or A\or B\or C\or D\or E\or F\or G\or H\or I\or J\or
   K\or L\or M\or N\or O\or P\or Q\or R\or S\or T\or U\or V\or W\or X\or
   Y\or Z\else\@ctrerr\fi}%
  \renewcommand\theequation{\@Alph\c@section.\@arabic\c@equation}%
\section*{APPENDIX \@Alph\c@section}}


\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{0.4\p@}
\setlength\doublerulesep{0.4\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{0.4\p@}
\renewcommand\theequation{\@arabic\c@equation}
\chardef\@ne=1
\chardef\tw@=2
\chardef\thr@@=3
\def\hb@xt@{\hbox to}
\def\eqnarray{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip 0.5\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip 0.5\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr
}
\def\endeqnarray{%
      \@@eqncr
      \egroup
      \global\advance\c@equation\m@ne
   $$\global\@ignoretrue
}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
           {\@float{figure}}
               {\end@float}
               
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\noindent\fontsize{8}{10}\selectfont #1 #2}%
  \ifdim \wd\@tempboxa >\hsize
    \noindent\fontsize{8}{10}\bfseries\sffamily\selectfont #1\hspace*{7pt} #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}



\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\xdef\bio@figure{figure}
               

\newcounter{tblcap}
\renewcommand\thetblcap{\@arabic\c@tblcap}

\newsavebox{\tabbox}
\newsavebox{\tabnotesbox}
\newlength{\tablen}
\newlength{\tabnoteslen}
\newcommand\centerlast{%
  \advance\leftskip by 0pt plus 1fil%
  \advance\rightskip by 0pt plus -1fil%
   \parfillskip0pt plus 2fil\relax}



\long\def\tblcaption#1#2{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
%\addtolength{\tablen}{-18pt}
\vskip8pt\refstepcounter{tblcap} 
\noindent\parbox[c]{\textwidth}{
\fontsize{9}{11}\bfseries\sffamily
\noindent\mbox{\noindent\centerlast
 TABLE \thetblcap}\hskip4.5pt  {\fontsize{9}{11}\bfseries\sffamily #1}}\\[-6.75pt]
\fontsize{8}{10}\selectfont \mbox{#2}}



\long\def\tblcaptionnotes#1#2#3{
\sbox{\tabbox}{#2}
\sbox{\tabnotesbox}{\fontsize{8}{10}\selectfont #2}
\settowidth{\tablen}{\usebox{\tabbox}}
\settowidth{\tabnoteslen}{\usebox{\tabnotesbox}}
%\addtolength{\tablen}{-18pt}
\vskip8pt\refstepcounter{tblcap} 
\noindent\parbox[c]{\textwidth}{
\fontsize{8}{10}\sffamily\bfseries\selectfont 
\noindent\mbox{\noindent
 TABLE \thetblcap}\hskip8pt  {\fontsize{8}{10}\bfseries\sffamily #1}}\\[-8.7pt]
\fontsize{8}{10}\selectfont \mbox{#2}\\[5pt]
\parbox[@{}l]{20pc}{\fontsize{8}{10}\selectfont #3}
}



\long\def\lantblcaption#1#2{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
%height=\the\textheight
%\addtolength{\tablen}{-18pt}
\vskip8pt\refstepcounter{tblcap} 
\centering\parbox[c]{48pc}{
%\fontsize{9}{11}\selectfont 
\centerlast{\noindent\mbox{\noindent\centerlast
 Table \thetblcap.}\hskip4.5pt  {\it #1}}}\\[7.75pt]
\fontsize{9}{11}\selectfont \mbox{#2}}


\long\def\lantblcaptionnotes#1#2#3{
\sbox{\tabbox}{#2}
\sbox{\tabnotesbox}{\fontsize{9}{11}\selectfont #2}
\settowidth{\tablen}{\usebox{\tabbox}}
\settowidth{\tabnoteslen}{\usebox{\tabnotesbox}}
\addtolength{\tabnoteslen}{-2.5pt}
%\addtolength{\tablen}{-18pt}
\vskip8pt\refstepcounter{tblcap} 
\centering\parbox[c]{48pc}{
%\fontsize{9}{11}\selectfont 
\centerlast\noindent\mbox{\noindent\centerlast
 Table \thetblcap.}\hskip4.5pt  {\it #1}}\\[7.75pt]
\fontsize{9}{11}\selectfont \mbox{#2}\\[5pt]
\hspace*{4pt}\parbox[@{}l]{\tabnoteslen}{\fontsize{8}{10}\selectfont #3}
}


\def\tblhead#1{\hline\\[-9pt]#1\\\hline\\[-9.75pt]}
\def\lastline{\\\hline}


\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\RequirePackage{latexsym}
\RequirePackage{amssymb}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\long\def\endstatement#1{\fontsize{9}{11}\selectfont\itshape\centering#1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD%%%%%%%%%%%%%%%%%%%%%%
% \newenvironment{thebibliography}[1]
%      {\section*{\refname\vspace*{-0.5pt}
%         %\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}
%         }%
%       \list{}%
%            {\labelwidth\z@%
%             \leftmargin 12\p@
%             \itemindent-\leftmargin
%             \@openbib@code
%             \usecounter{enumiv}%
%             \let\p@enumiv\@empty
%             \renewcommand\theenumiv{\@arabic\c@enumiv}}%
%             \fontsize{8}{9}\selectfont\parindent\z@
%             \itemsep5pt
%             \parskip0pt\relax%
%       \sloppy\clubpenalty4000\widowpenalty4000%
%       \sfcode`\.\@m
%           }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%%%%%%%%%%%%%%%%%%%%NEW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newdimen\bibindent
\setlength\bibindent{0em}%1.5em

\newenvironment{thebibliography}[1]
{\section*{\refname}\vspace*{3.5pt}%
      %\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \fontsize{8}{9}\selectfont\parindent\z@
            \advance\leftmargin \labelsep
            \advance\leftmargin by -4.5pt
            \labelsep 4pt %5
            \itemsep 2.7pt
            \parsep 0pt
            \listparindent 0pt%-10pt
            \itemindent \listparindent
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\. 1000}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist
      \endstatement 
     \vspace*{3pt}}
\newcommand\newblock{\hskip0.11em\@plus0.33em\@minus0.07em}
\let\@openbib@code\@empty
\renewcommand*{\@biblabel}[1]{#1.}%%%%%Ambika
% \renewcommand\footnoterule{%
%   \kern-2\p@
%   \hrule\@width.4\columnwidth
%   \kern2.6\p@}
  %%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%  
    
\newcommand\@makefntext[1]{%
    \parindent 6pt%
    \noindent
{\hb@xt@12pt{\hss\@makefnmark}
  \vspace*{-4pt} \fontsize{8}{8}\selectfont{#1}}}
\addtolength{\skip\footins}{6pt}    
\newcommand\refname{\fontsize{10}{12}\sffamily REFERENCES and FOOTNOTES}
\newcommand\figurename{FIGURE}
\newcommand\tablename{Table}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{ABSTRACT}
\newcommand\today{}
\edef\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\pagestyle{headings}%%%headings%%%Ambika
\pagenumbering{arabic}
\AtBeginDocument{\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\setcounter{footnote}{2}
}

% 
\def\addcontentsline#1#2#3{}%
%%\addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}}
\let\@period=. \def\period{\hbox{\@period}}
\mathcode`\.="8000
{\catcode`\.=\active\gdef.{\bio@period}}
\def\bio@period{\hbox{\@period}}






 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin_CROPMARKS

\newdimen\markwidth
\newdimen\marklength
\newdimen\markdistance
\newdimen\trimwidth
\newdimen\trimheight
\newdimen\rhsink
\newdimen\leftmarginodd
\newdimen\leftmargineven
\newdimen\leftmarginodda
\newdimen\leftmarginevena
\newdimen\auxaaa
\newdimen\auxbbb
\newdimen\auxccc
\newdimen\auxddd
\newdimen\tempbbb
\newdimen\rhsinka


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numbers supplied down here!
\markwidth=.5pt %0.4 truept % Thickness of the mark line
\marklength=20 truept % Length of the mark line
\markdistance=-6 truept % Distance of the real corner to the beg. of the mark
\trimwidth=202.22mm
\trimheight=276.2mm
\rhsink=4.55pc %%%48 truept % Sink from the top of paper to the BASE of RH
                 %  where the crops are hooked to (2pc + height of 10pt digit
\leftmarginodd=5pc %%%77pt %%%  yyyy
\leftmargineven=4.76pc %%77pt %%%

\leftmarginodda=5pc %%%77pt %%%  yyyy
\leftmarginevena=4.76pc %%77pt %%%
\rhsinka=5.1pc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\auxaaa\marklength\advance\auxaaa\markdistance
 \def\oddevenswitch{\ifodd\c@page \auxbbb=\leftmarginodd
  \else\auxbbb=\leftmargineven\fi}%

\def\oddevenswitcha{\ifodd\c@page \tempbbb=\leftmarginodda
  \else\tempbbb=\leftmarginevena\fi}%


\def\Vtrim{\VRHDW{20pt}{0pt}{1pt}}%
\def\Htrim{\VRHDW{1pt}{0pt}{20pt}}%

\def\ULtrim{\koo{-\markdistance}{0pt}{\zeroLR{\Htrim}}%
  \koo{0pt}{-\markdistance}{\zeroLR{\Vtrim}}}%
\def\URtrim{\koo{\markdistance}{0pt}{\zeroLL{\Htrim}}%
  \koo{0pt}{-\markdistance}{\zeroLL{\Vtrim}}}%
\def\LLtrim{\koo{-\markdistance}{0pt}{\zeroUR{\Htrim}}%
  \koo{0pt}{\markdistance}{\zeroUR{\Vtrim}}}%
\def\LRtrim{\koo{\markdistance}{0pt}{\zeroUL{\Htrim}}%
  \koo{0pt}{\markdistance}{\zeroUL{\Vtrim}}}%

\def\CROPMARKS{\oddevenswitch
  \koo{-\auxbbb}{-\rhsink}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}

\def\CROPMARKSA{\oddevenswitcha
  \koo{-\tempbbb}{-\rhsinka}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}



\def\CROPMARKS{\oddevenswitch
  \koo{-\auxbbb}{-\rhsink}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}

\def\CCROPMARKS{\oddevenswitch
 \hbox to0pt{\kern-\auxbbb\vbox to0bp{\kern-\rhsink\hbox to\trimwidth{%
 \koo{-\auxaaa}{0pt}{\VRHDW{\markwidth}{0pt}{\marklength}}%
 \koo{-\markwidth}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}%
 \hfill
 \koo{\markdistance}{0pt}{\VRHDW{\markwidth}{0pt}{\marklength}}%
 \koo{0pt}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}}%
 \nointerlineskip\vskip\trimheight
 \nointerlineskip\hbox to\trimwidth{%
 \koo{-\auxaaa}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
 \koo{-\markwidth}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}%
 \hfill
 \koo{\markdistance}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
 \koo{0pt}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}}\vss}\hss}}%

%\def\CROPMARKS{} %%%%%%%%% no cropmarks

%%%%%%%%%b-KOO
\long\def\koo#1#2#3{\vbox to0pt{\hsize0pt\kern #2\hbox to0pt{\kern
#1{#3}\hss}\vss}}
\long\def\zeroCC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}}}
\long\def\zeroUC#1{\vbox to0pt{\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}}}
\long\def\zeroCR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}\vss}}
\long\def\zeroUR#1{\vbox to0pt{\hbox to0pt{\hss #1}\vss}}
\long\def\zeroLL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}}}
\long\def\zeroCL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}\vss}}
\long\def\zeroUL#1{\vbox to0pt{\hbox to0pt{#1\hss}\vss}}
\def\VRHDW#1#2#3{\vrule height #1 depth #2 width #3}
\def\HRHDW#1#2#3{\hrule height #1 depth #2 width #3}
\def\HD#1#2{\vrule height #1pt depth #2pt width 0pt\relax}
%%%%%%%%%e-KOO


\def\cropson#1{
\let\cropmarkonoff=#1
\ifx N\cropmarkonoff
\def\CROPMARKSA{}
\def\CROPMARKS{}
\fi
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end_CROPMARKS

%\pagestyle{headings}

%%%%%%%%%%%%%%%%%%%%Theorem Definition

 \newtheoremstyle{theorem}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
 \theoremstyle{theorem}
\newtheorem{theorem}{\sc Theorem}[section]
%\renewcommand\thetheorem{\arabic{theorem}}

 \newtheoremstyle{lemma}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
 \theoremstyle{lemma}
 \newtheorem{lemma}[theorem]{\sc Lemma}
 

\newtheoremstyle{example}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
\theoremstyle{example}
\newtheorem{example}[theorem]{\sc Example}


\newtheoremstyle{Assumption}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
\theoremstyle{Assumption}
\newtheorem{Assumption}[theorem]{\sc Assumption}

\newtheoremstyle{ Proposition}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
\theoremstyle{ Proposition}
\newtheorem{Proposition}{\sc Proposition}

\newtheoremstyle{remark}{12pt}{12pt}{\rm}{}{\sffamily}{ }{ }{}
\theoremstyle{remark}
\newtheorem{remark}[theorem]{\sc Remark}

\renewenvironment{proof}[1][Proof]{\noindent\textit{#1. } }{\hfill$\square$}

\def\smlc{\fontencoding{OT1}\fontfamily{ptm}\fontseries{m}\fontshape{sc}\fontsize{10pt}{12pt}\selectfont}
% \newenvironment{references}
%                { \vbox{\vspace*{24pt}\centering{\smlc References}}\vspace*{-13pt}
%                           \fontsize{9}{11}\selectfont\list{}{  
%                          \topsep0pt
%                          \listparindent -12pt%
%                          \leftmargin 12pt
%                          \itemindent    \listparindent
%                          \itemsep 4.5pt
%                         \parsep        \z@ \@plus\p@}%
%                 \item\relax}
%                {\endlist\vspace*{11pt}}



\newenvironment{references}
     {\section*{\refname\vspace*{-0.5pt}
        %\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}
        }%
      \list{}%
           {\labelwidth\z@%
            \leftmargin 2\p@
            \itemindent-\leftmargin
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
            \fontsize{8}{9}\selectfont\parindent\z@
            \itemsep5pt
            \parskip0pt\relax%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m
          }
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist
     }

\sfcode`\. 1000





\def\cropmark#1{
\let\crop=#1
\ifx N\crop
\def\CROPMARKSA{}
\def\CROPMARKS{}
\fi
}





\def\op{{\rm [}}
\def\cp{{\rm ]}}

  \let\geq=\geqslant
  \let\leq=\leqslant
  \let\ge=\geqslant
  \let\le=\leqslant
\jno{kxi000}

%%%%%%%%%%%%%updated 4thAug-- Ambika
%\newcommand\doi[1]{doi: #1}
\def\doi#1{\gdef\@doinum{#1}}
%%2.7pc

\def\doiline{\vbox to 3pc{\vskip 15pt
\hfill{\vrule width 2.165in height 0.4pt}
 
\vskip 0pt \hfill {\fontsize{7}{9}\selectfont\copyright}\ {\fontsize{8}{10}\selectfont
2006 by the Biophysical Society}\\[-14pt]

\hfill {\fontsize{8}{10}\selectfont\@doinum}}}

\def\addrdagger{{\fontsize{5.65}{6}$^{\mbox{\textdagger}}$}}
\def\authdagger{{\fontsize{7.07}{8}$^{\mbox{\textdagger}}$}}
\frenchspacing

