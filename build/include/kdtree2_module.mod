  +2  r   k820309    s          18.0        4�`                                                                                                          
       /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/src/enkf/kdtree2.f90 KDTREE2_MODULE              KDTREE2 KDTREE2_RESULT TREE_NODE KDTREE2_CREATE KDTREE2_DESTROY KDTREE2_N_NEAREST KDTREE2_N_NEAREST_AROUND_POINT KDTREE2_R_NEAREST KDTREE2_R_NEAREST_AROUND_POINT KDTREE2_SORT_RESULTS KDTREE2_R_COUNT KDTREE2_R_COUNT_AROUND_POINT KDTREE2_N_NEAREST_BRUTE_FORCE KDTREE2_R_NEAREST_BRUTE_FORCE                                                     
                            @                              
                         @                               '                    #DIS    #IDX                 � $                                             	                � $                                                                  @                              'P                    #HEAP_SIZE    #ELEMS                � $                                                                                                                              0                �$                                                              #KDTREE2_RESULT              &                                                                                        	                                                         &         @                                
     P                      #RESULTS_IN    #PQ                                                                                 &                                           #KDTREE2_RESULT    %         @                                                   	       #A    #DIS    #IDX                                                   P               #PQ              
                                      	                
                                             %         @                                                   	       #A    #DIS    #IDX                                                   P               #PQ              
                                      	                
                                                               @                               '                    #DIMEN    #N    #THE_DATA    #IND    #SORT    #REARRANGE    #REARRANGED_DATA    #ROOT                � $                                                                                                                              0                � $                                                                                                                             0               �$                                                         	            &                   &                                                                                 y	                                                          �$                                          h                            &                                                                                 y                                                           � $                                   �                                                                                                             � $                                   �                                                                                                            �$                                         �                	            &                   &                                                                                 y	                                                           �$                                   p                   #TREE_NODE                                          y#TREE_NODE                                                                  �  @                              'p              	      #CUT_DIM    #CUT_VAL    #CUT_VAL_LEFT     #CUT_VAL_RIGHT !   #L "   #U #   #LEFT $   #RIGHT %   #BOX &                � D                                                              � D                                            	                � D                                             	                � D                             !               	                � D                              "                               � D                              #                               �D                              $     p                     #TREE_NODE                 �D                              %     p                      #TREE_NODE               �D                              &            (              	     #INTERVAL '             &                                                                                 y#INTERVAL '                                                                @  @                          '     '                    #LOWER (   #UPPER )                �                              (                	                �                              )               	   &         @                                *                           #INPUT_DATA +   #DIM ,   #SORT -   #REARRANGE .   #KDTREE2              `                              +                   	 	              &                   &                                                     
 @                               ,                     
 @                               -                     
 @                               .           #         @                                   /                    #TP 0             D P                               0                    #KDTREE2    #         @                                   1                    #TP 2   #QV 3   #NN 4   #RESULTS 5             D P                               2                    #KDTREE2              
                                 3                   	              &                                                     
  @                               4                     D `                               5                                   &                                           #KDTREE2_RESULT    #         @                                   6                    #TP 7   #IDXIN 8   #CORRELTIME 9   #NN :   #RESULTS ;             D P                               7                    #KDTREE2              
                                  8                     
                                  9                     
  @                               :                     D `                               ;                                   &                                           #KDTREE2_RESULT    #         @                                   <                    #TP =   #QV >   #R2 ?   #NFOUND @   #NALLOC A   #RESULTS B             D P                               =                    #KDTREE2              
                                 >                   	              &                                                     
                                 ?     	                D @                               @                      
  @                               A                     D `                               B                                   &                                           #KDTREE2_RESULT    #         @                                   C                    #TP D   #IDXIN E   #CORRELTIME F   #R2 G   #NFOUND H   #NALLOC I   #RESULTS J             D P                               D                    #KDTREE2              
                                  E                     
                                  F                     
                                 G     	                D @                               H                      
  @                               I                     D `                               J                                   &                                           #KDTREE2_RESULT    #         @                                  K                    #NFOUND L   #RESULTS M             
  @                               L                     D `                               M                   '                &                                           #KDTREE2_RESULT    %         @                                 N                           #TP O   #QV P   #R2 Q             D P                               O                    #KDTREE2              
                                 P                   	              &                                                     
                                 Q     	      %         @                                 R                           #TP S   #IDXIN T   #CORRELTIME U   #R2 V             D P                               S                    #KDTREE2              
                                  T                     
                                  U                     
                                 V     	      #         @                                   W                    #TP X   #QV Y   #NN Z   #RESULTS [             D P                               X                    #KDTREE2              
@ @                              Y                   	 !             &                                                     
                                  Z                     D                                 [                   "                &                                           #KDTREE2_RESULT    #         @                                   \                    #TP ]   #QV ^   #R2 _   #NFOUND `   #RESULTS a             D P                               ]                    #KDTREE2              
@ @                              ^                   	 $             &                                                     
                                 _     	                D @                               `                      D@                               a                   %                &                                           #KDTREE2_RESULT                  �  @                          b     '�                   #DIMEN c   #NN d   #NFOUND e   #BALLSIZE f   #CENTERIDX g   #CORRELTIME h   #NALLOC i   #REARRANGE j   #OVERFLOW k   #QV l   #RESULTS m   #PQ n   #DATA o   #IND p                � D                              c                                � D                              d                               � D                              e                               � D                             f               	               � D                              g                                                                               �              999                � D                              h                                                                               '              9999                 � D                              i                               � D                              j                               � D                              k             	                  �D                             l            (              
   	            &                                                      �D                              m            p                    #KDTREE2_RESULT              &                                                        � D                              n     P       �              #PQ                �D                             o                            	            &                   &                                                       �D                              p            h                            &                                              �   \      fn#fn $   �   0  b   uapp(KDTREE2_MODULE )   ,  @   J   KDTREE2_PRECISION_MODULE .   l  @   J   KDTREE2_PRIORITY_QUEUE_MODULE =   �  b       KDTREE2_RESULT+KDTREE2_PRIORITY_QUEUE_MODULE A     H   a   KDTREE2_RESULT%DIS+KDTREE2_PRIORITY_QUEUE_MODULE A   V  H   a   KDTREE2_RESULT%IDX+KDTREE2_PRIORITY_QUEUE_MODULE 1   �  j       PQ+KDTREE2_PRIORITY_QUEUE_MODULE ;     �   a   PQ%HEAP_SIZE+KDTREE2_PRIORITY_QUEUE_MODULE 7   �  �   a   PQ%ELEMS+KDTREE2_PRIORITY_QUEUE_MODULE    U  p       R_SINGLE+KINDS 8   �  h       PQ_CREATE+KDTREE2_PRIORITY_QUEUE_MODULE C   -  �   a   PQ_CREATE%RESULTS_IN+KDTREE2_PRIORITY_QUEUE_MODULE 8   �  i       PQ_INSERT+KDTREE2_PRIORITY_QUEUE_MODULE :   6  P   a   PQ_INSERT%A+KDTREE2_PRIORITY_QUEUE_MODULE <   �  @   a   PQ_INSERT%DIS+KDTREE2_PRIORITY_QUEUE_MODULE <   �  @   a   PQ_INSERT%IDX+KDTREE2_PRIORITY_QUEUE_MODULE =     i       PQ_REPLACE_MAX+KDTREE2_PRIORITY_QUEUE_MODULE ?   o  P   a   PQ_REPLACE_MAX%A+KDTREE2_PRIORITY_QUEUE_MODULE A   �  @   a   PQ_REPLACE_MAX%DIS+KDTREE2_PRIORITY_QUEUE_MODULE A   �  @   a   PQ_REPLACE_MAX%IDX+KDTREE2_PRIORITY_QUEUE_MODULE    ?	  �       KDTREE2    �	  �   a   KDTREE2%DIMEN    �
  �   a   KDTREE2%N !   :    a   KDTREE2%THE_DATA    F  �   a   KDTREE2%IND    :  �   a   KDTREE2%SORT "   �  �   a   KDTREE2%REARRANGE (   �    a   KDTREE2%REARRANGED_DATA    �  �   a   KDTREE2%ROOT    \  �       TREE_NODE "     H   !   TREE_NODE%CUT_DIM "   _  H   !   TREE_NODE%CUT_VAL '   �  H   !   TREE_NODE%CUT_VAL_LEFT (   �  H   !   TREE_NODE%CUT_VAL_RIGHT    7  H   !   TREE_NODE%L      H   !   TREE_NODE%U    �  _   !   TREE_NODE%LEFT     &  _   !   TREE_NODE%RIGHT    �    !   TREE_NODE%BOX    �  f       INTERVAL    �  H   a   INTERVAL%LOWER    C  H   a   INTERVAL%UPPER    �  �       KDTREE2_CREATE *     �   a   KDTREE2_CREATE%INPUT_DATA #   �  @   a   KDTREE2_CREATE%DIM $   �  @   a   KDTREE2_CREATE%SORT )   >  @   a   KDTREE2_CREATE%REARRANGE     ~  P       KDTREE2_DESTROY #   �  U   a   KDTREE2_DESTROY%TP "   #  m       KDTREE2_N_NEAREST %   �  U   a   KDTREE2_N_NEAREST%TP %   �  �   a   KDTREE2_N_NEAREST%QV %   q  @   a   KDTREE2_N_NEAREST%NN *   �  �   a   KDTREE2_N_NEAREST%RESULTS /   Q  �       KDTREE2_N_NEAREST_AROUND_POINT 2   �  U   a   KDTREE2_N_NEAREST_AROUND_POINT%TP 5   &  @   a   KDTREE2_N_NEAREST_AROUND_POINT%IDXIN :   f  @   a   KDTREE2_N_NEAREST_AROUND_POINT%CORRELTIME 2   �  @   a   KDTREE2_N_NEAREST_AROUND_POINT%NN 7   �  �   a   KDTREE2_N_NEAREST_AROUND_POINT%RESULTS "   �  �       KDTREE2_R_NEAREST %     U   a   KDTREE2_R_NEAREST%TP %   `  �   a   KDTREE2_R_NEAREST%QV %   �  @   a   KDTREE2_R_NEAREST%R2 )   ,  @   a   KDTREE2_R_NEAREST%NFOUND )   l  @   a   KDTREE2_R_NEAREST%NALLOC *   �  �   a   KDTREE2_R_NEAREST%RESULTS /   L  �       KDTREE2_R_NEAREST_AROUND_POINT 2   �  U   a   KDTREE2_R_NEAREST_AROUND_POINT%TP 5   9   @   a   KDTREE2_R_NEAREST_AROUND_POINT%IDXIN :   y   @   a   KDTREE2_R_NEAREST_AROUND_POINT%CORRELTIME 2   �   @   a   KDTREE2_R_NEAREST_AROUND_POINT%R2 6   �   @   a   KDTREE2_R_NEAREST_AROUND_POINT%NFOUND 6   9!  @   a   KDTREE2_R_NEAREST_AROUND_POINT%NALLOC 7   y!  �   a   KDTREE2_R_NEAREST_AROUND_POINT%RESULTS %   "  a       KDTREE2_SORT_RESULTS ,   z"  @   a   KDTREE2_SORT_RESULTS%NFOUND -   �"  �   a   KDTREE2_SORT_RESULTS%RESULTS     Z#  h       KDTREE2_R_COUNT #   �#  U   a   KDTREE2_R_COUNT%TP #   $  �   a   KDTREE2_R_COUNT%QV #   �$  @   a   KDTREE2_R_COUNT%R2 -   �$  {       KDTREE2_R_COUNT_AROUND_POINT 0   ^%  U   a   KDTREE2_R_COUNT_AROUND_POINT%TP 3   �%  @   a   KDTREE2_R_COUNT_AROUND_POINT%IDXIN 8   �%  @   a   KDTREE2_R_COUNT_AROUND_POINT%CORRELTIME 0   3&  @   a   KDTREE2_R_COUNT_AROUND_POINT%R2 .   s&  m       KDTREE2_N_NEAREST_BRUTE_FORCE 1   �&  U   a   KDTREE2_N_NEAREST_BRUTE_FORCE%TP 1   5'  �   a   KDTREE2_N_NEAREST_BRUTE_FORCE%QV 1   �'  @   a   KDTREE2_N_NEAREST_BRUTE_FORCE%NN 6   (  �   a   KDTREE2_N_NEAREST_BRUTE_FORCE%RESULTS .   �(  y       KDTREE2_R_NEAREST_BRUTE_FORCE 1   )  U   a   KDTREE2_R_NEAREST_BRUTE_FORCE%TP 1   o)  �   a   KDTREE2_R_NEAREST_BRUTE_FORCE%QV 1   �)  @   a   KDTREE2_R_NEAREST_BRUTE_FORCE%R2 5   ;*  @   a   KDTREE2_R_NEAREST_BRUTE_FORCE%NFOUND 6   {*  �   a   KDTREE2_R_NEAREST_BRUTE_FORCE%RESULTS #   +  �       TREE_SEARCH_RECORD )   ,  H   !   TREE_SEARCH_RECORD%DIMEN &   X,  H   !   TREE_SEARCH_RECORD%NN *   �,  H   !   TREE_SEARCH_RECORD%NFOUND ,   �,  H   !   TREE_SEARCH_RECORD%BALLSIZE -   0-  �   !   TREE_SEARCH_RECORD%CENTERIDX .   �-  �   !   TREE_SEARCH_RECORD%CORRELTIME *   .  H   !   TREE_SEARCH_RECORD%NALLOC -   �.  H   !   TREE_SEARCH_RECORD%REARRANGE ,   /  H   !   TREE_SEARCH_RECORD%OVERFLOW &   W/  �   !   TREE_SEARCH_RECORD%QV +   �/  �   !   TREE_SEARCH_RECORD%RESULTS &   �0  X   !   TREE_SEARCH_RECORD%PQ (   �0  �   !   TREE_SEARCH_RECORD%DATA '   �1  �   !   TREE_SEARCH_RECORD%IND 