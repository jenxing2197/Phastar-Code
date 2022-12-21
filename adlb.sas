/******************************************************************************
 *
 *Program Name      : adlb.sas
 *Protocol/Study    : Deciphera – MOTION –DCC3014 
 *Type              : ADaM
 *Description       : Program to create the ADLB (Laboratory Test Results, Analysis) Dataset
 *
 *Author            : Douglas Saberio
 *Date created      : 20Jun2022
 *Input datasets    : 
 *Macro used        : %checklog
 *Notes             : 
 *
 ******************************************************************************
 * Amended By       : Josh Enxing   
 * Date Amended     : 13DEC2022
 * Amendment        : Added code to filter by visits in SV according to DMC request
 ******************************************************************************
 * Amended By       : Josh Enxing   
 * Date Amended     : 14DEC2022
 * Amendment        : Fixing issues with multiple records per subject/visit/param with anl01fl=Y, incorrectly categorizing some visits, including visit not in SV
 ******************************************************************************
 * Amended By       : Josh Enxing   
 * Date Amended     : 16DEC2022
 * Amendment        : Fixing issues with erroneously including records with missing values
 ******************************************************************************
 * Amended By       : Josh Enxing   
 * Date Amended     : 19DEC2022
 * Amendment        : Updating anl01fl derivation to only be 'Y' for latest non-missing post-baseline value for a visit when there is a baseline reading
 ******************************************************************************/

/*Clear work data sets, work and output sessions*/
proc datasets lib = work nolist memtype = data kill;
run;

quit;

dm output 'clear';
dm log 'clear';

****************************************************************************
 formats for AVISITS, This can be modfified later
****************************************************************************;
data safetyfmta;
   length label $50 start 8;
   retain fmtname 'AVIS';
   start=-1;
   label='Screening';
   end=start;

   output;
   start=1;
   label='Cycle 1 Day 1';
   end=start;

   output;
   start=1.15;
   label='Cycle 1 Day 15';
   end=start;

   output;

   do start=4 to 8;
      vis1=start-2;
      label=catx('','Cycle',cat(vis1),' Day 1');
      end=vis1;

      output;
   end;

   start=6.99;
   label='End of Part 1';
   end=start;

   output;

   do start=10 to 16;
      vis1=start-3;
      label=catx('','Cycle',cat(vis1),' Day 1');
      end=vis1;

      output;
   end;

   start=25;
   label='Cycle 1 Day 1/Cycle 1 Day 1 Crossover';
   end=start;

   output;
   start=26;
   label='Cycle 1 Day 15/Cycle 1 Day 15 Crossover';
   end=start;

   output;

   do start=27 to 37;
      vis1=start-24;
      label=catx('','Cycle',cat(vis1),' Day 1','/Cycle',cat(vis1,','),' Day 1 Crossover');
      end=start;

      output;
   end;

   start=44;

   do label="Follow-Up","Follow-Up Week 12","Follow-Up Week 28","Follow-Up Week 44","Follow-Up Week 60"
      ,"Follow-Up Week 76","Follow-Up Week 92","Follow-Up Week 108";
      start+1;
      end=start;

      output;
   end;

   start=99;
   end=99;

   label= 'End of Treatment';
   output;
run;

data safetyfmt;
   length label $40 start 8;
   retain fmtname 'AVIST';
   start=-1;
   label='Screening';
   end=start;

   output;
   start=1;
   label='C1D1';
   end=start;

   output;
   start=1.15;
   label='C1D15';
   end=start;

   output;

   do start=4 to 8;
      vis1=start-2;
      label=catt('C',cat(vis1),'D1');
      end=vis1;

      output;
   end;

   start=6.99;
   label='EOP1';
   end=start;

   output;

   do start=10 to 16;
      vis1=start-3;
      label=catt('C',cat(vis1),'D1');
      end=vis1;

      output;
   end;

   start=25;
   label='C1D1/C1D1CO';
   end=start;

   output;
   start=26;
   label='C1D15/C1D15CO';
   end=start;

   output;

   do start=27 to 32;
      vis1=start-24;
      label=catx('','C',cat(vis1,','),'D1','/C',cat(vis1,','),'D1CO');
      end=start;

      output;
   end;

   start=44;

   do label="FUP","FUPW12","FUPWk28","FUP Wk44","FUP Wk60","FUP Wk76","FUP Wk92","FUP Wk108";
      start+1;
      end=start;

      output;
   end;

   start=99;
   end=99;

   label= 'EOT';
   output;
run;

data safetyfmt1(rename=(start1=start label1=label fmtname1=fmtname));
   length fmtname1 start1 $40 label1 8;
   set safetyfmt;
   start1=label;

   *end1=start1;
   label1=end;
   fmtname1='$AVISN';
   drop start end label fmtname vis1;
run;

data safetyfmt;
   set safetyfmt;
   start=end;
run;

data safetyfmta;
   set safetyfmta;
   start=end;
run;

proc format library=work cntlin=safetyfmt;
run;

proc format library=work cntlin=safetyfmt1;
quit;

proc format library=work cntlin=safetyfmta;
run;

****************************************************************************
---------------- End of Formats ------------------------------
****************************************************************************;

************************* ADSL **************************;
data _adsl;
   set adam.adsl;
   keep studyid usubjid subjid siteid tr: cutoffdt scrnfl saffl ittfl randdt ENRLFL region tumorloc dctdt1 prrther prsther digdy dislocl 
      dislocs;
run;

************************* SV data **************************;
proc sort data = raw.sv out = sv nodupkey;
   by subject instancename;
run;

data visits;
   set sv;
   subjid = strip(subject);
   visit = strip(instancename);

   if VISYN_STD = 'Y';
   keep subjid visit;
run;

************************* Macro to read in the different LB datasets and derive LBCAT for merging **************************;
%macro dsin(din,var)/ minoperator;

   data p_&din.(rename=( studyid=_studyid siteid=_siteid visit=visitx lbcat=lbcatx) drop= project: enviro: TARGETDAYS DATAPAGEID PAGEREPEATNUMBER recordid recordposition MINCREATED MAXUPDATED savets);
      set raw.lb_&din.;
      lbcat=upcase(scan(datapagename,findw(upcase(datapagename), "&var", ' ', 'E')));
      subjid=subject;
      visit=instancename;
      visitnumx=int(folderseq);

      %if %upcase(&din.) eq LBU and %upcase(&var.) eq URINALYSIS %then
         %do;
            array varx{*} HGB_LBORRES GLU_LBORRES KET_LBORRES PRO_LBORRES SPG_LBORRES_RAW;

            if lbperf = 'YES';

            do i=1 to dim(varx);
               lbtest=vlabel(varx(i));
               lbtestcd=scan(vname(varx(i)),1,'_');
               lborres=varx(i);

               if i<5 and lborres = '' then
                  lborres = 'MISSING';

               if ^missing(lborres) then
                  output;
            end;

            drop HGB_: GLU_: KET_: PRO_: SPG_:;
         %end;
      %else %if %upcase(&din.) ne (PG) and %upcase(&din.) ne (LBUM) %then
         %do;
            if upcase(lbperf) eq 'YES';
         %end;
   run;

   *JE-13DEC2022: added merge in since prod was not filtering by visits in SV;
   proc sort data = p_&din.;
      by subjid visitx;
   run;

   data _lb_&din.;
      merge visits (in=a rename=(visit=visitx)) p_&din. (in=b);
      by subjid visitx;

      if a and b;
   run;

%mend dsin;

%dsin(chem,CHEMISTRY);
%dsin(coag,COAGULATION);
%dsin(hema,HEMATOLOGY);
%dsin(lbu,URINALYSIS);
%dsin(lbum,URINALYSIS);
%dsin(pg,PREGNANCY);
%dsin(ser,SEROLOGY);

*------------- Append all datasets -------------------;
data all_set;
   length visit  lbcat $40;
   set _lb_: indsname=x;
   dsnam=x;

   if dsnam eq 'WORK._LB_LBUM' then
      do;
         if index(upcase(lbtest),'CAST')>0 or index(upcase(lbtest),'LEUKOCYTES') > 0 then
            lborres = strip(put(lborres1,best.));
         else lborres=strip(coalescec(lbumoro,lbumor));

         if lbumoro ne '' then
            lborresc = strip(lbumoro);
         else lborresc = strip(lbumor);

         if index(lbtest,'OTHER')>0 then
            lbtest = lbtestoth;
      end;
   else lborres=strip(coalescec(lborres,lborres1_raw));
   lbcat=lbcatx;
   visit=visitx;

   if lbdat = . and recorddate ne . then
      do;
         put 'WAR' 'NING: recorddate replacing lbdat for ' visitx ' for ' subjid;
         lbdat = recorddate;
      end;

   keep visit: lbcat: subject subjid lborres: lbtestcd lbtest lbdat lbdat_: lbtim foldername folderseq datapagename lbperf dsnam pg:;
run;

proc sort data= all_set;
   by subjid lbcat visit;
run;

*LK-26OCT2022: Update to use dcc_301403001_lab_prod_16sep2022 instead of the lb _ext copied from DMC1 dataset;
proc sort data= raw.dcc_301403001_lab_prod_16sep2022 out=dcc_301403001_lab_prod_16sep2022(drop=siteid);
   *nodupkey dupout=dup;
   by subjid visit;
run;

proc sql;
   alter table visits
      modify subjid char(20) format=$20.
      modify visit char(40) format=$40.;
quit;

data dcc;
   merge dcc_301403001_lab_prod_16sep2022 (in=a) visits (in=b);
   by subjid visit;

   if a and b;
run;

proc sort data = dcc;
   by subjid lbcat visit accsnnum lbtestcd;
run;

************************* Merge in Central data with Local data **************************;
data all_lb(rename=subjid_=subjid);
   length lbcat $150 subjid_ $8;
   merge all_set dcc;
   by subjid lbcat visit;

   if upcase(lbcat)='PREGNANCY' then
      do;
         p_lborres=lborres;
         lbdat = .;
      end;

   subjid_=subjid;
   drop subjid;

   *********************************************************
     Remove HCG information that has been collected in long format
   • For HCG, we are getting 2 result lines per visit – 1 being “positive/negative” and one being “numeric” (<0.6 for all). 
   Do we need to keep both of these, and if so, which one is needed for the baseline? (if any, as it’s only done at the screening visit.)
   **************************************************************************;
   if lbtestcd eq 'HCG' and ^missing(RPTU) then
      delete;
run;

data lb;
   attrib    
      STUDYID  label = 'Study Identifier'                         length = $21   
      USUBJID  label = 'Unique Subject Identifier'                length = $30   
      SUBJID   label = 'Subject Identifier for the Study'         length = $8    
      SITEID   label = 'Study Site Identifier'                    length = $5
      avalc length = $40
      anrind length = $40 
      lbtype length=$20 
      LBORRESU label = 'Original Units'                           length = $40
      LBORRES  label = ''                           length = $200
      PARCAT1  label = 'Parameter Category'                       length = $20   
      PARAMCD  label = 'Parameter Code'                           length = $200  
   ;
   set all_lb(rename=lborresu=lborresu_);
   by subjid lbcat visit;
   studyid= "DCC-3014-03-001";
   usubjid=catx("-",studyid,subjid);
   siteid=scan(subjid,1,"-");
   parcat1=lbcat;
   lbtox= toxgr;

   *******URINALYSIS **********;
   if ^missing(lbdat) then
      adtx=datepart(lbdat);

   if ^missing(lbdat) and ^missing(lbtim) then
      do;
         if length(lbtim) eq 4 then
            lbtim=cat(0,lbtim);
         v=catx('',put(adtx,yymmdd10.),lbtim);
         adtmx=input(v, e8601dt.);

         if ^missing(adtmx) then
            do;
               atm=timepart(adtmx);
               adt=datepart(adtmx);
            end;
      end;

   *******PREGNANCY **********;
   if ^missing(pgdat) then
      adtx=datepart(pgdat);

   if ^missing(pgdat) and ^missing(pgtim) then
      do;
         if strip(pgtim) ne '' then
            do;
               v=catx('',put(adtx,yymmdd10.),pgtim);
               adtmx=input(v,e8601dt.);
            end;

         if ^missing(adtmx) then
            do;
               atm=timepart(adtmx);
               adt=datepart(adtmx);
            end;
      end;

   *******OTHERS **********;
   if length(lbdtm)=10 then
      adty=input(lbdtm,yymmdd10.);
   else if length(lbdtm)=19 then
      do;
         adtmy=input(lbdtm,e8601dt.);
         atm=timepart(adtmy);
         adty=datepart(adtmy);
      end;

   if upcase(parcat1)='URINALYSIS' then
      do;
         avisit=visit;
         avisitn=visitnumx;

         if dsnam eq 'WORK._LB_LBUM' then
            do;
               if index(upcase(lbtest),'CAST') > 0 or index(upcase(lbtest),'LEUKOCYTES') > 0 then
                  avalc = upcase(lborres);
               else avalc = lborresc;
            end;
         else avalc = upcase(lborres);
         lbtype='LOCAL';

         ******-------  Convert character to numeric ---------**********;
         if length(avalc) eq length(compress(avalc,'.','kd')) and missing(aval) then
            aval=input(avalc,best. -L);
         param=lbtest;

         if upcase(param) eq 'SPECIFY GRAVITY (CHARACTER)' then
            param='Specific Gravity';
         else if upcase(param) eq 'KETONE' then
            param='Ketones';
         else if upcase(param) eq 'HYAL CAST' then
            param='Hyaline Casts';
         else if upcase(param) eq 'RBC CASTS' then
            param='RBC Casts';
         else if upcase(param) eq 'WBC CASTS' then
            param='WBC Casts';
         else param=propcase(param);

         if upcase(param) = 'SPECIFIC GRAVITY' then
            lbtestcd = 'SPGRAV';
         else if upcase(param) = 'GLUCOSE' then
            lbtestcd = 'GLUC';
         else if upcase(param) = 'LEUKOCYTES' then
            lbtestcd = 'WBC';
         else if upcase(param) = 'PROTEIN' then
            lbtestcd = 'PROT';
         else if upcase(param) = 'KETONES' then
            lbtestcd = 'KETONES';
         else if upcase(param) = 'HEMOGLOBIN' then
            lbtestcd = 'HGB';
         else if upcase(param) = 'BACTERIA' then
            lbtestcd = 'BACT';
         else if upcase(param) = 'HYALINE CASTS' then
            lbtestcd = 'CSHYAL';
         else if upcase(param) = 'RBC CASTS' then
            lbtestcd = 'CSRBC';
         else if upcase(param) = 'WBC CASTS' then
            lbtestcd = 'CSWBC';
         else if upcase(param) = 'MUCOUS' then
            lbtestcd = 'MUC';
         else if upcase(param) = 'SQUAMOUS EPITHELIAL CELLS' then
            lbtestcd = 'EPISQCE';
         else if upcase(param) = 'TRANSITIONAL EPITHELIAL CELLS' then
            lbtestcd = 'EPITCE';
         paramcd='U'||lbtestcd;
         adt=coalesce(adtx,adty);
         adtm=coalesce(adtmx,adtmy);
      end;
   else if upcase(parcat1)='PREGNANCY' then
      do;
         avisit=visit;
         avisitn=visitnumx;
         avalc=upcase(p_lborres);
         parcat1=lbcat;
         paramcd='PG';
         param='Pregnancy Test Result';
         adt=adtx;
         adtm=adtmx;
         lbtype='LOCAL';
      end;
   else
      do;
         avisit=visit;
         avisitn=input(visitnum,best.);
         avalc=upcase(siresc);
         parcat1=lbcat;
         paramcd=lbtestcd;
         anrhi=input(sinrhi,best.);
         anrlo=input(sinrlo,best.);
         a1hi=input(sinrhi,best.);
         a1lo=input(sinrlo,best.);
         lbspid=accsnnum;
         lborres=rptresc;
         lborresu=rptu;
         lbstresc=siresc;
         lbstresn=siresn;
         lbstresu=siu;
         lbornrlo=rptnrlo;
         lbornrhi=rptnrhi;
         lbstnrlo=sinrlo;
         lbstnrhi=sinrhi;
         adt=coalesce(adty,adtx);
         adtm=coalesce(adtmy,adtmx);
         lbtype='CENTRAL';

         *** Converting LBORRES to numeric form;
         * If lborres begins with (>, <, >=, <=, =>, =<) then only take numeric portion of result;
         * Elsewise, if full portion of result is numeric then take full result;
         if siresn^=. then
            aval=siresn;
         else if index(avalc,'<')>0 or index(avalc,'>')>0 then
            do;
               aval=input(compress(avalc,'<>=','s'),best.);
            end;

         if siu^='' then
            param=catx(' (',lbtest,strip(siu)||")");
         else param=lbtest;

         if param="Bilirubin (umol/L)" and upcase(parcat1)="CHEMISTRY" THEN
            param="Bilirubin, Total (umol/L)";

         * anrind;
         if not missing(lbstresc) then
            do;
               if substr(lbstresc,1,1) not in ('<' '>') then
                  do;
                     if ^cmiss(lbstnrhi,lbstnrlo) and (input(lbstnrlo,best.) <= lbstresn <= input(lbstnrhi,best.)) then
                        anrind = "NORMAL";
                     else if ^missing(lbstnrlo) and not missing(lbstresn) and lbstresn < input(lbstnrlo,best.) then
                        anrind = "LOW";
                     else if ^missing(lbstnrhi) and not missing(lbstresn) and lbstresn > input(lbstnrhi,best.) then
                        anrind = "HIGH";
                  end;
               else if substr(lbstresc,1,1) = '<' /*and index(substr(lbstresc,2),'+')=0*/
                  then

                     do;
                        if input(substr(lbstresc,2),??best.) <= input(lbstnrlo,best.) then
                           anrind="LOW";
                        else if input(lbstnrhi,best.) > input(substr(lbstresc,2),??best.) > input(lbstnrlo,best.) then
                           anrind="NORMAL";
                     end;
                  else if substr(lbstresc,1,1) = '>' /*and index(substr(lbstresc,2),'+')=0*/
                     then

                        do;
                           if input(substr(lbstresc,2),??best.) >= input(lbstnrhi,best.) then
                              anrind="HIGH";
                           else if input(lbstnrhi,best.) > input(substr(lbstresc,2),best.) > input(lbstnrlo,best.) then
                              anrind="NORMAL";
                        end;
            end;
      end;

   *JE-16DEC2022: added to remove rows with no data;
   if avalc ne '';
   format adt yymmdd10. atm time5. adtm e8601dt.;
run;

data lb1a x(keep=studyid usubjid param aval tr01sdt tr02sdt adt ady visit);
   length avisit $40 a1ind $10.;
   merge lb(in=a) _adsl(in=b);

   if a;
   by usubjid;

   ** Derive A1IND  **;
   if ~nmiss(aval,a1lo) & aval<a1lo then
      do;
         a1ind = 'LOW';

         if a1lo>0 then
            m_a1lo = aval/a1lo;
      end;
   else if ~nmiss(aval,a1hi) & aval>a1hi then
      do;
         a1ind = 'HIGH';

         if a1hi>0 then
            m_a1hi = aval/a1hi;
      end;
   else if ~nmiss(aval,a1lo, a1hi) & a1lo<=aval<=a1hi then
      a1ind = 'NORMAL';

   if ^missing(anrind) then
      a1ind=anrind;

   ** Derive R2A1HI **;
   if ^cmiss(a1hi,aval) and parcat1 eq 'CHEMISTRY' and lbtestcd in('ALT' 'AST' 'BILI') then
      r2a1hi=aval/a1hi;

   *
   r2a1hi=strip(put(r2a1hin,best.));

   ** Derive CRIT1A and CRIT1B and CRIT1FL **;
   length crit1 $60 crit1fl $1;

   * CRIT1;
   if lbtestcd in('ALT' 'AST') & ~nmiss(aval,a1hi) & m_a1hi>=3 then
      do;
         c1_1 = 1;
         crit1 = 'ALT or AST > 3XULN';
         crit1fl = 'Y';
      end;

   if lbtestcd='BILI' & ~nmiss(aval,a1hi) & m_a1hi>=2 then
      do;
         c1_2 = 1;
         crit1 = 'TBL > 2XULN';
         crit1fl = 'Y';
      end;

   if lbtestcd in('ALT' 'AST') then
      do;
         if 5< m_a1hi then
            do;
               crit2 = 'ALT or AST > 5XULN';
               critf2 ='Y';
            end;

         if 8< m_a1hi then
            do;
               crit3 = 'ALT or AST > 8XULN';
               critf3 ='Y';
            end;
      end;

   if ^cmiss(adt,tr01sdt) then
      ady=adt-tr01sdt+(adt>=tr01sdt);

   if (((tr01sdt<=adt and ^cmiss(adt,tr01sdt)) or cmiss(adt,tr01sdt) ne 2) or (tr01sdt<=adt<=tr02sdt and ^missing(tr02sdt))) then
      do;
         aperiod='Part 1';
         aperiodn=1;
      end;

   trtpn=trt01pn;
   trtp=trt01p;
   trta=trt01a;
   trtan=trt01an;

   if ^cmiss(adt,tr02sdt) and adt >=tr02sdt then
      do;
         trtpn=trt02pn;
         trtp=trt02p;
         trta=trt02a;
         trtan=trt02an;
         aperiod='Part 2';
         aperiodn=2;
      end;

   if adt<=tr01sdt and ^cmiss(adt,tr01sdt) then
      bf1fl='Y';

   if trt01a = 'Placebo' and trt02a = 'Vimseltinib' and adt<=tr02sdt and ^cmiss(adt,tr02sdt) then
      bf2fl='Y';
run;

proc sort data= lb1a(keep= usubjid lbtestcd adt c1_2 where=(c1_2 = 1)) out=c1_2 nodupkey;
   by usubjid lbtestcd adt;
run;

proc sort data= lb1a(keep= usubjid c1_1 lbtestcd adt where=(c1_1 = 1)) out=c1_1 nodupkey;
   by usubjid lbtestcd adt;
run;

proc sort data= lb1a;
   by usubjid lbtestcd adt;
run;

data c1_2;
   merge c1_2(in=a) c1_1(in=b);
   by usubjid lbtestcd adt;

   if a and b then
      anl03fl='Y';
run;

data lb1;
   merge lb1a(in=a) c1_2(in=b);
   by usubjid lbtestcd adt;

   if aperiod = 'Part 2' and trt01a = 'Placebo' and trt02a = 'Vimseltinib' then
      basetype = 'Part 2';
   else basetype = 'Part 1';
run;

proc sort data=lb1;
   by usubjid parcat1 bf1fl param adt;
run;

****************************
          two types of baseline one where lbcat in ('CHEMISTRY' 'HEMATOLOGY') and we have a numeric value and 
          where lbcat in ('CHEMISTRY' 'HEMATOLOGY') and we don't a have a numeric value
***************************;

*----------------- Baseline A -----------------------------;
data lb2a;
   set lb1(where=(lbcat in ('CHEMISTRY' 'HEMATOLOGY') and ^missing(aval)));
   by usubjid parcat1 bf1fl param adt;
   keep usubjid parcat1 param: bf1fl bf2fl aval: adt anrhi anrind a1ind tr01sdt tr02sdt;
run;

*----------------- Baseline B-----------------------------;
data lb2b;
   set lb1(where=(lbcat in ('CHEMISTRY' 'HEMATOLOGY') and ^missing(avalc)));
   by usubjid parcat1 bf1fl param adt;
   rename aval=avl avalc=avc anrhi=anrh a1ind=a1id;
   keep usubjid parcat1 param: bf1fl bf2fl aval: adt anrhi a1ind tr01sdt tr02sdt;
run;

data lb2c;
   merge lb2b(in=b) lb2a(in=a);
   by usubjid parcat1 bf1fl param adt;

   if (a and b) and ^missing(aval) then
      flag='Y';
   else if a and not b then
      flag1='Y';
   else flag2 ='Y';
run;

proc sort data= lb2c;
   by usubjid parcat1 bf1fl param flag flag1 flag2 adt;
run;

data lb2d;
   set lb2c;
   by usubjid parcat1 bf1fl param flag flag1 flag2 adt;

   if last.param and ^first.param then
      flag4='Y';
run;

*----------------- Baseline C-----------------------------;
data lb2e;
   set lb1(where=(lbcat ^in ('CHEMISTRY' 'HEMATOLOGY') and ^missing(avalc)));
   by usubjid parcat1 bf1fl param adt;

   if last.param then
      flag4='Y';
   keep usubjid parcat1 param: bf1fl bf2fl aval: adt anrhi anrind a1ind tr01sdt tr02sdt;
run;

data lb2f;
   set lb2d lb2e;
   by usubjid parcat1 bf1fl param;
   avalc=coalescec(avalc,avc);
   aval=coalesce(aval,avl);
   anrhi=coalesce(anrhi,anrh);
   a1ind=coalescec(a1ind,a1id);
run;

proc sort data= lb2f;
   by usubjid parcat1 bf1fl param flag flag1 flag2 flag4 adt aval avalc;
run;

data lb2;
   set lb2f;
   by usubjid parcat1 bf1fl param flag flag1 flag2 flag4 adt aval avalc;
   length basetype $40;

   if bf1fl='Y' and ^missing(tr01sdt) and last.param then
      do;
         ablfl='Y';
         bl01fl='Y';
         basetype='Part 1';
      end;
   else if bf2fl='Y' and ^missing(tr02sdt) and last.param then
      do;
         ablfl='Y';
         bl02fl='Y';
         basetype='Part 2';
         bl01fl='';
      end;

   if basetype = '' then
      basetype = 'Part 1';

   ******Add Urinalysis abnormal info if need be ********;
run;

data bl;
   length aperiod $6;
   set lb2;
   base=aval;
   basec=avalc;
   bnrind=anrind;
   b1ind=a1ind;
   aperiod = basetype;

   if ablfl='Y';
   keep usubjid parcat1 param: base basec bnrind bf1fl bf2fl basetype b1ind ablfl adt aperiod aval:;
run;

proc sort data=bl;
   by usubjid parcat1 param aperiod adt aval avalc;
run;

proc sort data=lb1;
   by usubjid parcat1 param aperiod adt aval avalc;
run;

data lb3a;
   merge lb1 (in=a) bl(in=b keep = usubjid parcat1 param aperiod adt ablfl aval avalc);
   by usubjid parcat1 param aperiod adt aval avalc;
run;

proc sort data = lb3a out = lb3b;
   by usubjid parcat1 param basetype;
run;

proc sort data = bl out = bla;
   by usubjid parcat1 param basetype;
run;

data lb3;
   merge bla(drop = ablfl) lb3b;
   by usubjid parcat1 param basetype;
run;

proc sort data=lb3;
   by usubjid parcat1 param adt;
run;

data lb4;
   set lb3;

   if ady>1 then
      do;
         if base^=. and aval^=. then
            chg=aval-base;

         if base not in (0 .) and aval^=. then
            pchg=round(((aval-base)/base)*100,0.01);
      end;

   /*lbrange*/
   if nmiss(anrlo,anrhi)=0 then
      lbrange = strip(put(anrlo,best.))||' - '||strip(put(anrhi,best.));
   else if missing(anrlo) and not missing(anrhi) then
      lbrange = '<'||strip(put(anrhi,best.));
   else if not missing(anrlo) and missing(anrhi) then
      lbrange = '>'||strip(put(anrlo,best.));
run;

proc sort data=lb4;
   by usubjid parcat1 param avisitn adt;
run;

data lb5;
   set lb4;
   by usubjid parcat1 param avisitn adt;

   if ^cmiss(adt,trtstdt) and ady>0 and ablfl ne 'Y' and last.avisitn and (base ne . or basec ne '') then
      anl01fl='Y';

   *
   keep trtstdt usubjid parcat1 param avisitn adt anl01fl ablfl;
run;

proc sort data=lb5;
   by studyid usubjid parcat1 lbspec paramcd adt adtm;
run;

data lb6(rename=(lbstnrhi_=lbstnrhi lbstnrlo_=lbstnrlo));
   set lb5(where=(^missing(param)));
   by studyid usubjid parcat1 lbspec paramcd adt adtm;
   lbstnrhi_ = input(lbstnrhi,best.);
   lbstnrlo_ = input(lbstnrlo,best.);

   *avisit = put(avisitn,avist.);
   if ady<1 then
      avisit="Screening";
   drop lbstnrlo lbstnrhi;
run;

*************************************************************************************************
      ------------------------ Derive AVISIT/N     ------------------------------------
*************************************************************************************************;

*Bring in Subject visit datasets;
data sv;
   attrib  
      STUDYID  label = 'Study Identifier'          length = $21 
      USUBJID  label = 'Unique Subject Identifier' length = $30;
   set raw.sv (rename=(studyid=_studyid))
      raw.sv2(rename=(studyid=_studyid))
      raw.sv3(rename=(studyid=_studyid))
      raw.sv6(rename=(studyid=_studyid))
      raw.sv7(rename=(studyid=_studyid))
      raw.sv8(rename=(studyid=_studyid));
   studyid="DCC-3014-03-001";
   usubjid=catx("-",studyid,subject);
   visit=strip(folder);
   visitnum=folderseq;

   if visdat^=. then
      svstdt=datepart(visdat);

   if upcase(visyn) = "YES" or upcase(svyn) = "YES";

   *if svstdt^=.;
   format svstdt date9.;
   keep studyid usubjid visit visitnum svstdt visyn svyn foldername;
run;

*Sort to keep one record per subject per visit;
proc sort data=sv;
   by studyid usubjid visitnum visit descending svstdt;
run;

data sva;
   length foldername $40;
   set sv(rename=(foldername=_foldername));
   by studyid usubjid visitnum visit descending svstdt;

   if first.visit;
   foldername=_foldername;
run;

data sv1;
   merge sva(in=a) _adsl(in=b keep=studyid usubjid randdt ENRLFL);
   by studyid usubjid;

   if b;

   /* Some subjects were randomized and did not start dose and are not in SV */
   if svstdt=. and visit in ("","SCREEN") then
      do;
         visit="Screening";
         visitnum=1;
         svstdt=randdt;
         foldername='Screening';
      end;

   output;

   if visit in ("C2D1","C3D1","C4D1","C5D1","C6D1","C7D1","C8D1","C9D1", "C10D1") then
      do;
         visit=strip(visit)||"3";
         visitnum=visitnum+0.5;

         if svstdt ne . then
            svstdt=svstdt+12;
         output;
      end;

   drop randdt;
run;

proc transpose data=sv1 out=t_sv(drop=_name_);
   by studyid usubjid;
   id visit;
   idlabel foldername;
   var svstdt;
run;

/**LK-Commented***/
data t_sva;
   set t_sv;

   if c1d15=. and c1d1^=. then
      c1d15=c1d1+14;

%macro loop;
   %do j= 2 %to 10;
      if c&j.d1=. and c%eval(&j.-1)d1^=. then
         c&j.d1=c%eval(&j.-1)d1+28;
   %end;
%mend loop;

%loop;
array visdt{*} c2d1 c3d1 c4d1 c5d1 c6d1 c7d1 c8d1 c9d1 c10d1;
array nvisdt{*} c2d13 c3d13 c4d13 c5d13 c6d13 c7d13 c8d13 c9d13 c10d13;

do i=1 to dim(visdt);
   if visdt{i}^=. then
      nvisdt{i}=visdt{i}+12;
end;

format c2d13 c3d13 c4d13 c5d13 c6d13 c7d13 c8d13 c9d13 c10d13 date9.;
drop i;
run;

proc transpose data=t_sva out=t_sv1a(rename=(col1=svstdt));
   by studyid usubjid;
   var screen c1d1 c1d15 c2d1 c2d13 c3d1 c3d13 c4d1 c4d13 c5d1 c5d13 c6d1 c6d13 c7d1 c7d13 c8d1 c8d13 c9d1 c9d13 c10d1 c10d13;
run;

data t_sv1;
   set t_sv1a;
   where svstdt^=.;
   length visit $50;
   visit=strip(_name_);
   foldername=_label_;
   drop _name_ _label_;
run;

proc sort data=lb6 out=visits1(keep=studyid usubjid adt dctdt1 ady visit visitnum avisit avisitn: tr02sdt trtstdt randdt) nodupkey;
   by studyid usubjid avisit adt;
run;

/*****/
data sv2;
   set sv;
run;

proc sort data=sva out = sv2_s;
   by studyid usubjid visit;
run;

proc sort data=t_sv1 out = t_sv1_s;
   by studyid usubjid visit foldername;
run;

data t_sv2;
   merge  t_sv1_s(in=a) sv2_s(in=b);
   by studyid usubjid visit foldername;

   *if b;
   avisit = visit;
   if visyn='YES';
run;

proc sort data = t_sv2;
  by usubjid svstdt;
run;

/**********/
proc sql noprint;
   create table visits11 as
      select a.*, b.svstdt, (b.svstdt-a.adt)+1 as diff,b.avisit as visitcd
         from visits1 as a right join t_sv2 as b
            on a.usubjid=b.usubjid and a.visit = b.foldername
         order by studyid, usubjid, adt, diff;
quit;

data visits11a;
   set visits11;
   by studyid usubjid adt diff;

   if visit ne 'Screening' and diff ne . then
      diff=abs(diff);

   if visit = 'Screening' and ady = 1 then 
      do;
	     avisit = 'Cycle 1 Day 1';
         visit = 'Cycle 1 Day 1';
         visitcd = 'C1D1';
         avisitn = 1.01; 
      end; 

   if visitcd in ("C2D13","C3D13","C4D13","C5D13","C6D13","C7D13","C8D13","C9D13") then
      do;
         visitcd=substr(visitcd,1, length(visitcd)-1);
      end;
run;

proc sort data=visits11a;
   by studyid usubjid adt diff;
run;

data visits11b;
   set visits11a;
   by studyid usubjid adt diff;

   if ady<1 then
      do;
         avisit="Screening";
         visitcd="Screening";
      end;

   drop diff svstdt ady;

   if visit eq 'End of Treatment' and (missing(tr02sdt) or missing(dctdt1)) then
      do;
         avisit='End of Treatment';
         visit=avisit;
      end;
   else if visit eq 'End of Treatment' and ^missing(tr02sdt) then
      do;
         visit='End of Treatment';
         ;
         visit=avisit;
      end;
   else if visit eq 'End of Part 1' then
      do;
         avisit='End of Part 1';
      end;

   if first.adt;
run;

proc sort data=lb6 out=lb6_s;
   by usubjid;
run;

data lb7;
   merge lb6_s(drop = studyid subjid siteid) _adsl;
   by usubjid;
run;

proc sort data=lb7 out=lb7_s;
   by studyid usubjid adt;
run;

data alllb(rename=(avisitx=avisit visitnum_x=visitnum));
   length avisit avisitcd avisitx $200 avisitn 8.;
   merge  lb7_s(in=b drop=  avisit avisitn) visits11b(in=a drop= visit visitnum:);
   by studyid usubjid adt;

   if a;

   if visit = 'Screening' then
      visitcd = 'Screening';
   test=put(visitcd,$AVISN.);
   avisitn=input(test,?? best.);
   avisitcd=visitcd;
   avisitx=avisit;
   visitx=put(avisitn,avis.);
   visitnum_x=input(visitnum,best.);

   if ablfl = 'Y' then
      do;
         avisitn=0;
         avisitx='Baseline';
         avisitcd = 'BASE';
      end;

   *keep studyid usubjid parcat1 lbspec paramcd adt adtm avis: vis:;
   drop avisit visitnum test;
run;

proc sort data = alllb;
  by usubjid parcat1 param avisitn adt aval avalc;
run;

data alllb1;
   set alllb;
   by usubjid parcat1 param avisitn adt aval avalc;

   if ^cmiss(adt,trtstdt) and ady>0 and ablfl ne 'Y' and last.avisitn and (base ne . or basec ne '') then
      anl01fl='Y';
run;


proc sort data=alllb1;
   by studyid usubjid parcat1 lbspec paramcd adt adtm aval avalc;
run;

data alllc;
   set alllb1 (where=(^missing(param)));
   by studyid usubjid parcat1 lbspec paramcd adt adtm aval avalc;

   if first.usubjid then
      aseq=0;
   aseq+1;

   *JE-13DEC2022: added to correct avisitn;
   if upcase(avisit) = 'SCREENING' then
      avisitn = -1;
   else if upcase(avisit) in ('EOP1', 'END OF PART 1') then
      avisitn = 6.99;
   else if upcase(avisit) in ('EOT', 'END OF TREATMENT') then
      avisitn = 99;
   else if upcase(avisit) = 'UNSCHEDULED' then
      avisitn = 98;
   else if index(upcase(avisit),'CYCLE')>0 and index(upcase(avisit),'DAY')=0 then
      avisitn = input(scan(strip(avisit),2),best.);
   else if index(upcase(avisit),'CYCLE')>0 and index(upcase(avisit),'DAY')>0 then
      avisitn = input(scan(strip(avisit),2),best.) + input(scan(strip(avisit),4),best.)/100;
   else if avisit ne 'Baseline' then
      put 'WAR' 'NING: visit name not coded for: ' avisit;
run;

*Template for attributes;
data lb_templ;
   attrib      
      STUDYID  label = 'Study Identifier'                         length = $21    
      USUBJID  label = 'Unique Subject Identifier'                length = $30    
      SUBJID   label = 'Subject Identifier for the Study'         length = $8 
      SITEID   label = 'Study Site Identifier'                    length = $5 
      LBSPID   label = 'Sponsor-Defined Identifier'               length = $40    
      LBSPEC   label = 'Specimen Type'                            length = $200   
      LBNAM    label = 'Vendor Name'                              length = $200   
      LBORRES  label = 'Result or Finding in Original Units'      length = $200   
      LBORRESU label = 'Original Units'                           length = $40    
      LBSTRESC label = 'Character Result/Finding in Std Format'   length = $200   
      LBSTRESN label = 'Numeric Result/Finding in Standard Units' length = 8  
      LBSTRESU label = 'Standard Units'                           length = $40    
      LBORNRLO label = 'Reference Range Lower Limit in Orig Unit' length = $40    
      LBORNRHI label = 'Reference Range Upper Limit in Orig Unit' length = $40    
      LBSTNRLO label = 'Reference Range Lower Limit-Std Units'    length = 8  
      LBSTNRHI label = 'Reference Range Upper Limit-Std Units'    length = 8  
      ASEQ     label = 'Sequence Number'                          length = 8  
      TRTP     label = 'Planned Treatment for Period'             length = $40    
      TRTPN    label = 'Planned Treatment for Period (N)'         length = 8  
      TRTA     label = 'Actual Treatment for Period'              length = $40    
      TRTAN    label = 'Actual Treatment for Period (N)'          length = 8  
      PARCAT1  label = 'Parameter Category'                       length = $20    
      PARAMCD  label = 'Parameter Code'                           length = $200   
      PARAM    label = 'Parameter'                                length = $200   
      AVAL     label = 'Analysis Value'                           length = 8  
      AVALC    label = 'Analysis Value (C)'                       length = $200   
      BASE     label = 'Analysis Baseline Value'                  length = 8  
      BASEC    label = 'Analysis Baseline Value (C)'              length = $200   
      CHG      label = 'Change from Baseline'                     length = 8  
      BASETYPE label = 'Baseline Type'                            length = $200 
      VISIT   label = 'Visit'                                     length = $200   
      AVISITN  label = 'Analysis Visit (N)'                       length = 8  
      AVISIT   label = 'Analysis Visit'                           length = $200   
      VISITNUM  label = 'Analysis Visit (N)'                       length = 8  
      APERIOD  label = 'Period'                                   length = $200   
      APERIODN label = 'Period (n)'                               length = 8  
      ADT      label = 'Analysis Date'                            length = 8  
      ADTM     label = 'Analysis Datetime'                        length = 8  
      ADY      label = 'Analysis Relative Day'                    length = 8  
      ABLFL    label = 'Baseline Record Flag'                     length = $200   
      A1IND    label = 'Analysis Range 1 Indicator'               length = $10    
      A1HI     label = 'Analysis Range 1 Upper Limit'             length = 8  
      A1LO     label = 'Analysis Range 1 Lower Limit'             length = 8  
      ANRLO    label = 'Analysis Normal Range Lower Limit'        length = 8  
      ANRHI    label = 'Analysis Normal Range Upper Limit'        length = 8  
      ANRIND   label = 'Analysis Reference Range Indicator'       length = $40    
      B1IND    label = 'Baseline Analysis Range 1 Indicator '     length = $10    
      BNRIND   label = 'Baseline Reference Range Indicator'       length = $40  
      R2A1HI   label = 'Ratio to Analysis Range 1 Upper Limit'    length = 8  
      CRIT1    label = 'Analysis Criterion 1'                     length = $60    
      CRIT1FL  label = 'Criterion 1 Evaluation Result Flag'       length = $1 
      CRIT2    label = 'Analysis Criterion 2'                     length = $100   
      CRIT2FL  label = 'Criterion 2 Evaluation Result Flag'       length = $2 
      CRIT3    label = 'Analysis Criterion 3'                     length = $100   
      CRIT3FL  label = 'Criterion 3 Evaluation Result Flag'       length = $2 
      ANL01FL  label = 'Analysis 01 Flag'                         length = $200   
      TRTSDT   label = 'Start Date of First Treatment'            length = 8  
      TRTSTDT  label = 'Date of First Exposure to Treatment'      length = 8  
      TRTENDT  label = 'Date of Last Exposure to Treatment'       length = 8  
      TR01SDT  label = 'Date of First Exposure in Period 01'      length = 8  
      TR01EDT  label = 'Date of Last Exposure in Period 01'       length = 8  
      TR02SDT  label = 'Date of First Exposure in Period 02'      length = 8  
      TR02EDT  label = 'Date of Last Exposure in Period 02'       length = 8  
      CUTOFFDT label = 'Data Cutoff Date'                         length = 8  
      SCRNFL   label = 'Screening Failure Flag'                           length = $2 
      SAFFL    label = 'Safety Population Flag'                   length = $2 
      ITTFL    label = 'Intent-to-Treat Population Flag'          length = $2 
      REGION   label = 'Geographic Region of Site'                length = $200    
      tumorloc label = 'Tumor Location'                           length = $200    
      PRRTHER  label = 'Prior TGCT radiation therapy'             length = $200   
      PRSTHER  label = 'Prior Systemic Therapy'                   length = $200   
      DIGDY    label = 'Time from diagnosis to first dose date'   length = 8  
      DISLOCL  label = 'Disease located in large joints '         length = $200   
      DISLOCS  label = 'Disease located in small joints '         length = $200       
      LBTOX    label = 'Toxicity'                                 length = $200   
      LBTYPE   label = 'Type of Laboratory'                       length = $20    
      ANL03FL  label = 'Analysis Flag 03'                         length = $1;
   call missing(of _all_);
   stop;
run;

data adam.adlb;
   retain studyid usubjid subjid siteid lbspid lbspec lbnam lbtype lborres lborresu lbstresc lbstresn lbstresu lbornrlo lbornrhi   
      lbstnrlo lbstnrhi aseq trtp trtpn trta trtan parcat1 paramcd param aval avalc base basec chg basetype /*visit visitnum*/
   avisitn avisit avisitcd aperiod aperiodn adt
   adtm ady ablfl a1ind a1hi a1lo anrlo anrhi anrind b1ind bnrind  r2a1hi  crit1 crit1fl crit2 crit2fl crit3 crit3fl anl01fl trtsdt trtstdt 
   trtendt tr01sdt tr01edt tr02sdt tr02edt cutoffdt scrnfl saffl ittfl region tumorloc prrther prsther digdy dislocl dislocs basetype lbtox
   anl03fl
   ;
   set lb_templ alllc;
   keep studyid usubjid subjid   siteid lbspid lbspec lbnam  lborres lborresu lbstresc lbstresn lbstresu lbornrlo lbornrhi   
      lbstnrlo lbstnrhi aseq trtp trtpn trta trtan parcat1 paramcd param aval avalc base basec chg basetype /*visit visitnum avisitcd*/
   avisitn avisit  aperiod aperiodn adt
   adtm ady ablfl a1ind a1hi a1lo anrlo anrhi anrind b1ind bnrind  r2a1hi  crit1 crit1fl crit2 crit2fl crit3 crit3fl anl01fl trtsdt trtstdt 
   trtendt tr01sdt tr01edt tr02sdt tr02edt cutoffdt scrnfl saffl ittfl region tumorloc prrther prsther digdy dislocl dislocs basetype lbtox
   lbtype anl03fl;
run;

/*proc sort data=sandbox.adlb_qc out=adlb_qc; by studyid usubjid parcat1 lbspec paramcd param adt adtm;run;*/
/*proc sort data=adam.adlb out=adlb; by studyid usubjid parcat1 lbspec paramcd param adt adtm;run;*/
/*proc compare base=adlb compare=adlb_qc listall; id studyid usubjid parcat1 lbspec paramcd param adt adtm;quit;*/
%checklog;
