  �B  �   k820309    s          18.0        T�`                                                                                                          
       /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/src/gsi/berror.f90 BERROR       ?       INIT_BERROR PCINFO CREATE_BERROR_VARS DESTROY_BERROR_VARS SET_PREDICTORS_VAR RESET_PREDICTORS_VAR INIT_RFTABLE INITABLE CREATE_BERROR_VARS_REG DESTROY_BERROR_VARS_REG QVAR3D NR NF VARPRD FPSPROJ BKGV_FLOWDEP FUT2PS NDX NDY NDX2 NMIX NYMX NFG NFNF NORM NXEM VPRECOND ADJUSTOZVAR DSSVS DSSV BKGV_WRITE BKGV_REWGTFCT HSWGT HZSCL BW PERT_BERR_FCT PERT_BERR NDEG NORH VS BL BL2 BE SLW2 SLW1 SLW MR INAXS WTXRS WTAXS NX NY INXRS JJ1 II2 JJ2 II JJ II1 TABLE ALV NHSCRF CWCOVEQQCOV                      @                              
       R_KIND I_KIND                      @                              
       NC3D NVARS MVARS                                                                                                                                                                                                                                                                                                                                                                                  #         @                                                                                                                          #         @                                   	                    #PCINFO%MAXSCAN 
                                                                                                                   
            #         @                                                       #CREATE_BERROR_VARS%LON1    #CREATE_BERROR_VARS%LAT1    #CREATE_BERROR_VARS%ITOTSUB    #CREATE_BERROR_VARS%IGLOBAL                                                                                                                                                                                                                                                                                 #         @                                                       #DESTROY_BERROR_VARS%NCLEN    #DESTROY_BERROR_VARS%NRCLEN                                                                                                              #         @                                                       #SET_PREDICTORS_VAR%NLON    #SET_PREDICTORS_VAR%NLAT    #SET_PREDICTORS_VAR%NSIG    #SET_PREDICTORS_VAR%LON2    #SET_PREDICTORS_VAR%LAT2    #SET_PREDICTORS_VAR%MAXSCAN    #SET_PREDICTORS_VAR%NCLEN    #SET_PREDICTORS_VAR%LON1    #SET_PREDICTORS_VAR%LAT1    #SET_PREDICTORS_VAR%ITOTSUB    #SET_PREDICTORS_VAR%IGLOBAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              #         @                                                       #RESET_PREDICTORS_VAR%NLON     #RESET_PREDICTORS_VAR%NLAT !   #RESET_PREDICTORS_VAR%NSIG "   #RESET_PREDICTORS_VAR%LON2 #   #RESET_PREDICTORS_VAR%LAT2 $   #RESET_PREDICTORS_VAR%MAXSCAN %   #RESET_PREDICTORS_VAR%NCLEN &   #RESET_PREDICTORS_VAR%LON1 '   #RESET_PREDICTORS_VAR%LAT1 (   #RESET_PREDICTORS_VAR%ITOTSUB )   #RESET_PREDICTORS_VAR%IGLOBAL *                                                                                                                                                                    !                                                     "                                                     #                                                     $                                                     %                                                     &                                                     '                                                     (                                                     )                                                     *            #         @                                   +     	              #INIT_RFTABLE%NLON ,   #INIT_RFTABLE%NLAT -   #INIT_RFTABLE%NSIG .   #INIT_RFTABLE%LON2 /   #INIT_RFTABLE%LAT2 0   #INIT_RFTABLE%LON1 1   #INIT_RFTABLE%LAT1 2   #INIT_RFTABLE%ITOTSUB 3   #INIT_RFTABLE%IGLOBAL 4   #MYPE 5   #RATE 6   #NNN 8   #SLI 9   #SLI1 <   #SLI2 >                                                                                         ,                                                     -                                                     .                                                     /                                                     0                                                     1                                                     2                                                     3                                                     4                      
                                 5                    
@ @                              6                    
    p          5 r 7       5 r 7                               
                                 8                    
  @                              9                    
        p        p        p         5 r :   5 r ;   p           5 r :   5 r ;     p          5 � p        r 8        5 r :   5 r ;     p          5 � p        r 8                              
 @                              <                    
        p        p        p         '   n                                           25 r =   n                                          1'   n                                          25 r =   n                                          1p           '   n                                      25 r =   n                                          1'   n                                          25 r =   n                                          1  p          5 � p        r 8        '   n                                      25 r =   n                                          1'   n                                          25 r =   n                                          1  p          5 � p        r 8                                                 
 @                              >                    
        p        p        p         '   n                                           25 r =   n                                          1'   n                                          25 r =   n                                          1p           '   n                                      25 r =   n                                          1'   n                                          25 r =   n                                          1  p          5 � p        r 8        '   n                                      25 r =   n                                          1'   n                                          25 r =   n                                          1  p          5 � p        r 8                                        #         @                                  ?                 
   #NXDIM @   #NYDIM A   #SLI B   #NTAX C   #IHWLB D   #IIX E   #JJX F   #FACTOR G   #TIN H   #IPOINT I                               
                                 @                     
                                 A                    
                                 B                    
 !       p        5 � p        r A   p        5 � p        r @   p          5 � p        r @     5 � p        r A     p            5 � p        r @     5 � p        r A     p                                    
  @                              C                     
                                 D                    D                                E                     #      p        5 � p        r @   p          5 � p        r @     5 � p        r A       5 � p        r @     5 � p        r A                              D                                F                     $      p        5 � p        r @   p          5 � p        r @     5 � p        r A       5 � p        r @     5 � p        r A                               
                                 G     
                
                                 H     
               
                                 I                     "   p          5 � p        r C       5 � p        r C                     #         @                                   J                    #CREATE_BERROR_VARS_REG%LON1 K   #CREATE_BERROR_VARS_REG%LAT1 L   #CREATE_BERROR_VARS_REG%ITOTSUB M   #CREATE_BERROR_VARS_REG%IGLOBAL N                                                                             K                                                     L                                                     M                                                     N            #         @                                   O                              @                                P                   
                &                   &                   &                                                      @                                Q                      @                                =                     @                                R                   
                &                                                      @                                S                       @                                 T                       @                                U                       @                                V                       @                                W                       @                                X                       @                                Y                       @                                Z                       @ @                              [                       @                                \                       @                                ]                       @                                ^                     @                                _                   
                &                                                      @                                 `                     @ @                              a                   
                &                   &                   &                                                    @ @                              b                   
                &                   &                   &                   &                                                      @                                 c                       @                                d     
                  @     �                           e                   
      p          & p        p            p                                     @ @   �                           f                   
      p          & p        p            p                                     @                                g     
                  @                                h     
                  @                                 i                      @ @                              7                       @                                j                       @                                k     
                @                                l                   
                &                                                    @                                m                   
                &                                                    @                                n                   
                &                                                    @                                o                   
                &                   &                                                    @                                p                   
                &                   &                                                    @                                q                   
                &                   &                                                      @                                r                     @                                s                                   &                   &                                                    @                                t                   
                &                   &                   &                                                    @                                u                   
                &                   &                   &                                                     @ @                              ;                      @ @                              :                     @                                v                                   &                   &                                                    @ @                              w                                   &                   &                   &                   &                                                    @ @                              x                                   &                   &                   &                   &                                                    @ @                              y                                   &                   &                   &                   &                                                    @ @                              z                                   &                   &                   &                   &                                                    @ @                              {                                   &                   &                   &                   &                                                    @ @                              |                                   &                   &                   &                   &                                                    @ @                              }                   
                &                   &                                                    @ @                              ~                   
                &                   &                   &                   &                                                      @                                                       @                                �               �   R      fn#fn    �   �  b   uapp(BERROR    �  N   J  KINDS     *  Q   J  CONTROL_VECTORS    {  p       R_KIND+KINDS    �  p       I_KIND+KINDS %   [  @       NC3D+CONTROL_VECTORS &   �  @       NVARS+CONTROL_VECTORS &   �  @       MVARS+CONTROL_VECTORS      �       INIT_BERROR    �  �       PCINFO /   H  @     PCINFO%MAXSCAN+RADINFO=MAXSCAN #   �  �       CREATE_BERROR_VARS 5   �  @     CREATE_BERROR_VARS%LON1+GRIDMOD=LON1 5   �  @     CREATE_BERROR_VARS%LAT1+GRIDMOD=LAT1 ;     @     CREATE_BERROR_VARS%ITOTSUB+GRIDMOD=ITOTSUB ;   C  @     CREATE_BERROR_VARS%IGLOBAL+GRIDMOD=IGLOBAL $   �  �       DESTROY_BERROR_VARS 6   
	  @     DESTROY_BERROR_VARS%NCLEN+JFUNC=NCLEN 8   J	  @     DESTROY_BERROR_VARS%NRCLEN+JFUNC=NRCLEN #   �	  $      SET_PREDICTORS_VAR 5   �  @     SET_PREDICTORS_VAR%NLON+GRIDMOD=NLON 5   �  @     SET_PREDICTORS_VAR%NLAT+GRIDMOD=NLAT 5   .  @     SET_PREDICTORS_VAR%NSIG+GRIDMOD=NSIG 5   n  @     SET_PREDICTORS_VAR%LON2+GRIDMOD=LON2 5   �  @     SET_PREDICTORS_VAR%LAT2+GRIDMOD=LAT2 ;   �  @     SET_PREDICTORS_VAR%MAXSCAN+RADINFO=MAXSCAN 5   .  @     SET_PREDICTORS_VAR%NCLEN+JFUNC=NCLEN 5   n  @     SET_PREDICTORS_VAR%LON1+GRIDMOD=LON1 5   �  @     SET_PREDICTORS_VAR%LAT1+GRIDMOD=LAT1 ;   �  @     SET_PREDICTORS_VAR%ITOTSUB+GRIDMOD=ITOTSUB ;   .  @     SET_PREDICTORS_VAR%IGLOBAL+GRIDMOD=IGLOBAL %   n  �      RESET_PREDICTORS_VAR 7   W  @     RESET_PREDICTORS_VAR%NLON+GRIDMOD=NLON 7   �  @     RESET_PREDICTORS_VAR%NLAT+GRIDMOD=NLAT 7   �  @     RESET_PREDICTORS_VAR%NSIG+GRIDMOD=NSIG 7     @     RESET_PREDICTORS_VAR%LON2+GRIDMOD=LON2 7   W  @     RESET_PREDICTORS_VAR%LAT2+GRIDMOD=LAT2 =   �  @     RESET_PREDICTORS_VAR%MAXSCAN+RADINFO=MAXSCAN 7   �  @     RESET_PREDICTORS_VAR%NCLEN+JFUNC=NCLEN 7     @     RESET_PREDICTORS_VAR%LON1+GRIDMOD=LON1 7   W  @     RESET_PREDICTORS_VAR%LAT1+GRIDMOD=LAT1 =   �  @     RESET_PREDICTORS_VAR%ITOTSUB+GRIDMOD=ITOTSUB =   �  @     RESET_PREDICTORS_VAR%IGLOBAL+GRIDMOD=IGLOBAL      �      INIT_RFTABLE /   �  @     INIT_RFTABLE%NLON+GRIDMOD=NLON /   �  @     INIT_RFTABLE%NLAT+GRIDMOD=NLAT /     @     INIT_RFTABLE%NSIG+GRIDMOD=NSIG /   [  @     INIT_RFTABLE%LON2+GRIDMOD=LON2 /   �  @     INIT_RFTABLE%LAT2+GRIDMOD=LAT2 /   �  @     INIT_RFTABLE%LON1+GRIDMOD=LON1 /     @     INIT_RFTABLE%LAT1+GRIDMOD=LAT1 5   [  @     INIT_RFTABLE%ITOTSUB+GRIDMOD=ITOTSUB 5   �  @     INIT_RFTABLE%IGLOBAL+GRIDMOD=IGLOBAL "   �  @   a   INIT_RFTABLE%MYPE "     �   a   INIT_RFTABLE%RATE !   �  @   a   INIT_RFTABLE%NNN !   �  d  a   INIT_RFTABLE%SLI "   S  (  a   INIT_RFTABLE%SLI1 "   {  (  a   INIT_RFTABLE%SLI2    �!  �       INITABLE    d"  @   a   INITABLE%NXDIM    �"  @   a   INITABLE%NYDIM    �"  t  a   INITABLE%SLI    X$  @   a   INITABLE%NTAX    �$  @   a   INITABLE%IHWLB    �$  $  a   INITABLE%IIX    �%  $  a   INITABLE%JJX      '  @   a   INITABLE%FACTOR    `'  @   a   INITABLE%TIN     �'  �   a   INITABLE%IPOINT '   T(  �       CREATE_BERROR_VARS_REG 9   G)  @     CREATE_BERROR_VARS_REG%LON1+GRIDMOD=LON1 9   �)  @     CREATE_BERROR_VARS_REG%LAT1+GRIDMOD=LAT1 ?   �)  @     CREATE_BERROR_VARS_REG%ITOTSUB+GRIDMOD=ITOTSUB ?   *  @     CREATE_BERROR_VARS_REG%IGLOBAL+GRIDMOD=IGLOBAL (   G*  H       DESTROY_BERROR_VARS_REG    �*  �       QVAR3D    K+  @       NR    �+  @       NF    �+  �       VARPRD    W,  @       FPSPROJ    �,  @       BKGV_FLOWDEP    �,  @       FUT2PS    -  @       NDX    W-  @       NDY    �-  @       NDX2    �-  @       NMIX    .  @       NYMX    W.  @       NFG    �.  @       NFNF    �.  @       NORM    /  @       NXEM    W/  �       VPRECOND    �/  @       ADJUSTOZVAR    #0  �       DSSVS    �0  �       DSSV    �1  @       BKGV_WRITE    �1  @       BKGV_REWGTFCT    32  �       HSWGT    �2  �       HZSCL    {3  @       BW    �3  @       PERT_BERR_FCT    �3  @       PERT_BERR    ;4  @       NDEG    {4  @       NORH    �4  @       VS    �4  �       BL    �5  �       BL2    6  �       BE    �6  �       SLW2    C7  �       SLW1    �7  �       SLW    �8  @       MR    �8  �       INAXS    o9  �       WTXRS    +:  �       WTAXS    �:  @       NX    ';  @       NY    g;  �       INXRS    <  �       JJ1    �<  �       II2    �=  �       JJ2    �>  �       II    [?  �       JJ    /@  �       II1    A  �       TABLE    �A  �       ALV    {B  @       NHSCRF    �B  @       CWCOVEQQCOV 