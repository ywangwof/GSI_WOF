  �  8   k820309    s          18.0        *�`                                                                                                          
       /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/src/gsi/lag_interp.f90 LAG_INTERP       	       LAG_ACCUR LAG_INIPGRID LAG_DELPGRID LAG_LOGCTE_P LAG_GRIDREL_IJK LAG_INDEX_H LAG_INT3D_NL LAG_INT3D_TL LAG_INT3D_AD                      @                              
       R_KIND I_KIND                      @                              
       NLON NLAT NSIG RLONS RLATS                      @                              
       ZERO ONE                  @                               '                                                                                                                                                                                                                                               @ @                                                     @ @                                                     @ @                              	                     @ @                              
                   
                &                                                    @ @                                                 
                &                                                                                            
                
                                 0.0                                                 
                
                       �?        1.0                                                                                                                                                                                                                                                                                                                                                                                    
       #         @                                                       #NEWGRID              
 @                                                 
              &                                           #         @                                                                @ @                                                 
                &                                           #         @                                                      #LON    #LAT    #P    #I    #J    #K              
                                      
                
                                      
                
  @                                   
                D @                                   
                 D @                                   
                 D @                                   
       %         @                                                           #LON !   #LAT "             
                                 !                     
                                 "           %         @                                #                    
       #FIELD $   #LON %   #LAT &   #P '   #LSPEC_I (   #LSPEC_R )                                    
                                 $                   
              &                   &                                                     
  @                              %     
                
  @                              &     
                
  @                              '     
                F @                              (                        p          p            p                                    F @                              )     
              
     p          p 
           p 
                         %         @                                *                    
       #LSPEC_I +   #LSPEC_R ,   #DLON -   #DLAT .   #DFIELD /             
  @                              +                       p          p            p                                    
                                 ,     
              
    p          p 
           p 
                                   
                                 -     
                
                                 .     
                
                                 /                   
              &                   &                                           #         @                                   0                    #LSPEC_I 1   #LSPEC_R 2   #ADINT3D 3   #ADLON 4   #ADLAT 5   #ADFIELD 6             
  @                              1                       p          p            p                                    
                                 2     
              
    p          p 
           p 
                                   
                                 3     
                
D                                4     
                 
D                                5     
                 
D                                6                   
               &                   &                                              �   Z      fn#fn     �   �   b   uapp(LAG_INTERP    ~  N   J  KINDS    �  [   J  GRIDMOD    '  I   J  CONSTANTS '   p  P       #UNLPOLY+ISO_C_BINDING    �  p       R_KIND+KINDS    0  p       I_KIND+KINDS    �  @       NLON+GRIDMOD    �  @       NLAT+GRIDMOD       @       NSIG+GRIDMOD    `  �       RLONS+GRIDMOD    �  �       RLATS+GRIDMOD    x  s       ZERO+CONSTANTS    �  s       ONE+CONSTANTS     ^  @       IGLOBAL+GRIDMOD     �  @       ITOTSUB+GRIDMOD    �  @       LAT2+GRIDMOD      @       LON2+GRIDMOD    ^  @       LAT1+GRIDMOD    �  @       LON1+GRIDMOD    �  @       LAG_ACCUR      U       LAG_INIPGRID %   s  �   a   LAG_INIPGRID%NEWGRID    �  H       LAG_DELPGRID    G	  �       LAG_LOGCTE_P     �	  v       LAG_GRIDREL_IJK $   I
  @   a   LAG_GRIDREL_IJK%LON $   �
  @   a   LAG_GRIDREL_IJK%LAT "   �
  @   a   LAG_GRIDREL_IJK%P "   	  @   a   LAG_GRIDREL_IJK%I "   I  @   a   LAG_GRIDREL_IJK%J "   �  @   a   LAG_GRIDREL_IJK%K    �  b       LAG_INDEX_H     +  @   a   LAG_INDEX_H%LON     k  @   a   LAG_INDEX_H%LAT    �  �       LAG_INT3D_NL #   P  �   a   LAG_INT3D_NL%FIELD !   �  @   a   LAG_INT3D_NL%LON !   4  @   a   LAG_INT3D_NL%LAT    t  @   a   LAG_INT3D_NL%P %   �  �   a   LAG_INT3D_NL%LSPEC_I %   H  �   a   LAG_INT3D_NL%LSPEC_R    �  �       LAG_INT3D_TL %   f  �   a   LAG_INT3D_TL%LSPEC_I %   �  �   a   LAG_INT3D_TL%LSPEC_R "   �  @   a   LAG_INT3D_TL%DLON "   �  @   a   LAG_INT3D_TL%DLAT $     �   a   LAG_INT3D_TL%DFIELD    �  �       LAG_INT3D_AD %   D  �   a   LAG_INT3D_AD%LSPEC_I %   �  �   a   LAG_INT3D_AD%LSPEC_R %   l  @   a   LAG_INT3D_AD%ADINT3D #   �  @   a   LAG_INT3D_AD%ADLON #   �  @   a   LAG_INT3D_AD%ADLAT %   ,  �   a   LAG_INT3D_AD%ADFIELD 