function [final_topology] = LPTD_Function_2(protein_data,helix_info,sheet_info) 



if helix_info==1
    
    
    pdb=pdbread(protein_data.name);  % online: PDBStruct = getpdb(PDBid)
    number_all_helix=numel(pdb.Helix);  % Number of Helices in the Constructed Model
    AtomNumber=numel(pdb.Model.Atom);  % Number of atoms in the Constructed Model
    start_residue_chain=[pdb.Helix.initSeqNum];  % Initial Sequence Number
    end_residue_chain=[pdb.Helix.endSeqNum];    % End Sequence Number
    length_helices=[pdb.Helix.length];   % length of helices
    resSeq=[pdb.Model.Atom.resSeq]; % Residues
    
    % chain
    k=1;
    helix_number=0;
    start_residue_chain=0;
    end_residue_chain=0;
    for i=1:number_all_helix
        if pdb.Helix(i).initChainID==protein_data.chain
            helix_number=helix_number+1;
            length_chain(k)=pdb.Helix(i).length;
            start_residue_chain(k)=pdb.Helix(i).initSeqNum;
            end_residue_chain(k)=pdb.Helix(i).endSeqNum;
            k=k+1;
        end
    end
    
    c=1;
    Matrix=zeros(helix_number,6);
    cc=1;
    for i=1:helix_number
        s=start_residue_chain(i);
        e=end_residue_chain(i);
        
        for j=1:AtomNumber
            sw=0;
            if ((pdb.Model.Atom(j).resSeq)==s) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                
                s_e_coordinate(1)=pdb.Model.Atom(j).X;
                s_e_coordinate(2)=pdb.Model.Atom(j).Y;
                s_e_coordinate(3)=pdb.Model.Atom(j).Z;
                
            end
            
            if ((pdb.Model.Atom(j).resSeq)==e) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                sw=1;
                s_e_coordinate(4)=pdb.Model.Atom(j).X;
                s_e_coordinate(5)=pdb.Model.Atom(j).Y;
                s_e_coordinate(6)=pdb.Model.Atom(j).Z;
                
                
            end
            if sw==1
                Matrix(cc,:)=s_e_coordinate; % Matrix for all start-end coordinate voxels
                cc=cc+1;
            end
        end
        
    end
    
    % Extract C_∝ voxeles of helices
    
    for i=1:helix_number
        m= Matrix(i,1);
        n= Matrix(i,2);
        o= Matrix(i,3);
        c=length_chain(i);
        k=1;
        l=0;  
        next=0;
        flag =1;
        while (flag)
            for j=1:AtomNumber
                if ((pdb.Model.Atom(j).X==m) && (pdb.Model.Atom(j).Y==n) && (pdb.Model.Atom(j).Z==o))
                    helix_matrix(k,1)=pdb.Model.Atom(j).X;
                    helix_matrix(k,2)=pdb.Model.Atom(j).Y;
                    helix_matrix(k,3)=pdb.Model.Atom(j).Z;
                    k=k+1;
                    next=start_residue_chain(i)+1;
                    l=l+1;
                end  % for first iteration
                
                if ((pdb.Model.Atom(j).resSeq)==next) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                    
                    helix_matrix(k,1)=pdb.Model.Atom(j).X;
                    helix_matrix(k,2)=pdb.Model.Atom(j).Y;
                    helix_matrix(k,3)=pdb.Model.Atom(j).Z;
                    k=k+1;
                    next=next+1;
                    l=l+1;
                    
                end
                
                if (l==c)
                    helix_voxels(i).Mat = helix_matrix;
                    helix_matrix = [];
                    flag = 0;
                    break;
                end % end if
            end % end for j=1:AtomNumber
        end % end while
    end % end for i=1:helix_number
    
    
    
    
    % Extract coordinate Voxels from  α-Sticks
    
    
    filename=protein_data.stick_hlces;
    d = csvread(filename);
    numberStick = unique (d(:,4)); % array of index
    cnt1 =1;
    n=size(d,1); % number of  rows
    stick_number_helices=size(numberStick);  % number of indexes
    stick_number_helices=stick_number_helices(1);
    
    for j=1:stick_number_helices
        value=numberStick(j); % select each #index
        k=1;
        for i =1:n
            if (d(i,4)==value)
                stick_Matrix(k,:)=d(i,1:3);
                k=k+1;
            end
        end
        stick_voxels_helices(j).Mat = stick_Matrix;
        stick_Matrix = [];
    end
    
    
    
    %%
    
    
    %**************************************************************************
    
    %     Construction of the complete weighted bipartite graph for α-helices  
    
    %**************************************************************************
    
    tic;
    
 %  Assign weights to the nodes of the graph by Improved Mahalanobis Squared Distance (IMSD)   
    
    for i = 1:helix_number
        k=1;
        for j =1:stick_number_helices
            A = helix_voxels(i).Mat;
            B = stick_voxels_helices(j).Mat;
            distance_mahal_AB = IMSD(A,B);
            min_distance = min(distance_mahal_AB);
            if min_distance>1e5
                min_distance=1e5
            end
            
            DD_AB(:,j) = min_distance;
            k=k+1;
        end
        d = [DD_AB];
        [minimum,index] = min(d(:));
        Mahalanobis(i).num_helix=i;
        Mahalanobis(i).num_stick=index;
        Mahalanobis(i).min=minimum;
        weight_matrix(i,:)=[d];   
    end
    
    
    
    % ********************* Linear programming Algorithm  *********************
  
    f=Lp_code(weight_matrix,0);  % Recall linear programing function
    
    for k=1:helix_number
        
        reduced_weight_matrix=weight_matrix;  % pre-assign stick=1 with all helices
        reduced_weight_matrix(:,1)=[];     %remove stick#1 (column=1)
        reduced_weight_matrix(k,:)=[]
        [lp_reduced,fval]=Lp_code(reduced_weight_matrix,0); % Run lp after remove stick#1 (column=1) and helix (row=k)
        fval_sum=fval+weight_matrix(k,1) ;
        
        % Adjust index of assignment by considering the removed column (column=1)
        
        for j=1:helix_number-1
            p=lp_reduced(j)
            if p==0         % index before remove row k
                lp_modify(j)=0
            end
            if p~=0
                lp_modify(j)=lp_reduced(j)+1
            end
        end
        
        w=1
        for i=1:helix_number
            if i~=k
                lp_topology(i).num_helix=i;
                lp_topology(i).num_stick=lp_modify(w);
                w=w+1;
            end
            
            if i==k
                lp_topology(i).num_helix=k;
                lp_topology(i).num_stick=1;
            end
        end
        
      
        
        lp_reduced_matrix(k).k=k;
        lp_reduced_matrix(k).lp_topology=[lp_topology];
        lp_reduced_matrix(k).lp_score=fval_sum;

        
    end
    
    % Sort topologies based on the minimum score of LP
    
    table_lp_matrix=struct2table(lp_reduced_matrix);
    sorted_table=sortrows(table_lp_matrix,'lp_score');
    sorted_table(:,1)=[];
    sorted_struct_lp=table2struct(sorted_table)
    
    % Assign rank for each topology
    for i=1:helix_number
        sorted_struct_lp(i).rank=i
        sorted_struct_lp(i).lp_topology=sorted_struct_lp(i).lp_topology;
        sorted_struct_lp(i).lp_score=sorted_struct_lp(i).lp_score;
    end
    
   
    
    % ******* Determine Top rank topology based on the lp_score to extract final assignment
    
    all_topologies=[sorted_struct_lp.lp_score]
    [minimum,rank_topology_helix] = min(all_topologies)
    final_topology_helix=sorted_struct_lp(rank_topology_helix).lp_topology
    
    
    
    % ********* Find the direction in topology based on the final assignment (DDA Algorithm)  *******************
    
    for i=1:helix_number
        
        if final_topology_helix(i).num_stick ~= 0    % If we have a matched pair (helix#, stick#)
            index=final_topology_helix(i).num_stick;
            len=length(stick_voxels_helices(index).Mat);  % Number of stick voxels
            c=1;
            inverse_stick_voxels = zeros(len,3);  % Preallocate for matrix
            for j=len:-1:1
                inverse_stick_voxels(c,:)=stick_voxels_helices(index).Mat(j,:);
                c=c+1;
            end
            
            A = helix_voxels(i).Mat;
            B1 = stick_voxels_helices(index).Mat;
            dtw_x1=dtw(A(:,1),B1(:,1));
            dtw_y1=dtw(A(:,2),B1(:,2));
            dtw_z1=dtw(A(:,3),B1(:,3));
            dtw_distance1 = sum(dtw_x1+dtw_y1+dtw_z1);
            
            B2 =inverse_stick_voxels;
            dtw_x2=dtw(A(:,1),B2(:,1));
            dtw_y2=dtw(A(:,2),B2(:,2));
            dtw_z2=dtw(A(:,3),B2(:,3));
            dtw_distance2 = sum(dtw_x2+dtw_y2+dtw_z2);
            
            
            if dtw_distance1<dtw_distance2 direction=1
            else direction=-1;
            end
                
        else        
            direction=0
            
        end
        final_topology_helix(i).Direction= direction;
       
        
    end
    
     final_topology=final_topology_helix;
     
end  % end helix info'

%%

% *************************************************************************************

% Extraction of geometrical features from generated model and cryo-EM map for sheets

% *************************************************************************************



if sheet_info==1
    
    pdb=pdbread(protein_data.name);  % online: PDBStruct = getpdb(PDBid)
    number_all_strand=numel(pdb.Sheet);  % Number of strands in in the Constructed Model
    AtomNumber=numel(pdb.Model.Atom);  % Number of atoms in the Constructed Model
    start_residue_chain=[pdb.Sheet.initSeqNum]; % Initial Sequence Number
    end_residue_chain=[pdb.Sheet.endSeqNum];    % End Sequence Number
    resSeq=[pdb.Model.Atom.resSeq]; % Residue 
    
    
    % Select chain`
    k=1;
    strand_number=0;
    start_residue_chain=0;
    end_residue_chain=0;
    for i=1:number_all_strand
        if pdb.Sheet(i).initChainID==protein_data.chain
            strand_number=strand_number+1;
            start_residue_chain(k)=pdb.Sheet(i).initSeqNum;
            end_residue_chain(k)=pdb.Sheet(i).endSeqNum;
            k=k+1;
        end
    end
    
    
    for i=1:strand_number
        length_strands(i) =end_residue_chain(i)-start_residue_chain(i)+1
        
    end
    
    
    c=1;
    Matrix=zeros(strand_number,6); 
    cc=1;
    for i=1:strand_number
        s=start_residue_chain(i);
        e=end_residue_chain(i);
        
        for j=1:AtomNumber
            sw=0;
            if ((pdb.Model.Atom(j).resSeq)==s) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                
                s_e_coordinate(1)=pdb.Model.Atom(j).X;
                s_e_coordinate(2)=pdb.Model.Atom(j).Y;
                s_e_coordinate(3)=pdb.Model.Atom(j).Z;
                
            end
            
            if ((pdb.Model.Atom(j).resSeq)==e) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                sw=1;
                s_e_coordinate(4)=pdb.Model.Atom(j).X;
                s_e_coordinate(5)=pdb.Model.Atom(j).Y;
                s_e_coordinate(6)=pdb.Model.Atom(j).Z;
                
                
            end
            if sw==1
                Matrix(cc,:)=s_e_coordinate; % Matrix start-end coordinate
                cc=cc+1;
            end
        end
        
    end
    
    for i=1:strand_number
        m= Matrix(i,1);
        n= Matrix(i,2);
        o= Matrix(i,3);
        c=length_strands(i);
        k=1;
        l=0;  % counter for the length of each strand
        next=0;
        flag =1;
        while (flag)
            for j=1:AtomNumber
                if ((pdb.Model.Atom(j).X==m) && (pdb.Model.Atom(j).Y==n) && (pdb.Model.Atom(j).Z==o))
                    strand_matrix(k,1)=pdb.Model.Atom(j).X;
                    strand_matrix(k,2)=pdb.Model.Atom(j).Y;
                    strand_matrix(k,3)=pdb.Model.Atom(j).Z;
                    k=k+1;
                    next=start_residue_chain(i)+1;
                    l=l+1;
                end  
                
                if ((pdb.Model.Atom(j).resSeq)==next) && (strcmp(pdb.Model.Atom(j).AtomName,'CA')==1)
                    
                    strand_matrix(k,1)=pdb.Model.Atom(j).X;
                    strand_matrix(k,2)=pdb.Model.Atom(j).Y;
                    strand_matrix(k,3)=pdb.Model.Atom(j).Z;
                    k=k+1;
                    next=next+1;
                    l=l+1;
                    
                end
                
                if (l==c)
                    strand_voxels(i).Mat = strand_matrix;
                    strand_matrix = [];
                    flag = 0;
                    break;
                end % End if
            end % End for j=1:AtomNumber
        end % End while 
    end % End for i=1:strand_number
    
  
    
    % **************** Generate More Voxels for strand Sticks *****************
    
    filename=protein_data.stick_strands;
    d = csvread(filename);
    numberStick = unique (d(:,4)); % array of index
    cnt1 =1;
    n=size(d,1); % number of all rows
    stick_number_strands=size(numberStick);  % number of indexes
    stick_number_strands=stick_number_strands(1);
    
    for j=1:stick_number_strands
        value=numberStick(j); % select each #index
        k=1;
        for i =1:n
            if (d(i,4)==value)
                stick_Matrix(k,:)=d(i,1:3);
                k=k+1;
            end
        end
        stick_voxels_strands(j).Mat = stick_Matrix;
        
        
        p1=stick_voxels_strands(j).Mat(1,:)
        p2=stick_voxels_strands(j).Mat(3,:)
        [x,y,z]= bresenham_line3d(p1,p2)
        
        generated_stick_voxels(j).Mat(:,1)=x
        generated_stick_voxels(j).Mat(:,2)=y
        generated_stick_voxels(j).Mat(:,3)=z
        
        stick_Matrix = [];
    end
    
    
    
    %%
    %**************************************************************************
    
    %     Construction of the weighted complete bipartite graph for sheets  
    
    %**************************************************************************
    
    tic;
    
 %  Assign weights to the nodes of the graph by Improved Mahalanobis Squared Distance (IMSD)   
    
    weight_matrix=[];
    for i = 1:strand_number
        k=1;
        for j =1:stick_number_strands
            A = strand_voxels(i).Mat;
            B = generated_stick_voxels(j).Mat;
            distance_mahal_AB = IMSD(A,B);
            min_distance = min(distance_mahal_AB);
            if min_distance>1e5
                min_distance=1e5
            end
            
            DD_AB_strand(:,j) = min_distance;
            k=k+1;
        end
        d_strand = [DD_AB_strand];
        [minimum,index] = min(d_strand(:));
        Mahalanobis(i).num_strand=i;
        Mahalanobis(i).num_stick=index;
        Mahalanobis(i).min=minimum;
        weight_matrix(i,:)=[d_strand];  
    end
    
    
    
    % ********************* Linear programming Algorithm  *********************
   
   
    f=Lp_code(weight_matrix,0);  % Recall linear programing function
    lp_topology=[];
    reduced_weight_matrix=[];
    lp_reduced=[];
    fval=0;
    lp_modify=[]
    lp_reduced_matrix=[]
    table_lp_matrix=[]
    sorted_struct_lp=[]
    all_topology_accuracy=[]
    
    for k=1:strand_number
        
        reduced_weight_matrix=weight_matrix; % pre-assign stick=1 with all helices
        reduced_weight_matrix(:,1)=[];     %remove stick#1 (column=1)
        reduced_weight_matrix(k,:)=[]
        [lp_reduced,fval]=Lp_code(reduced_weight_matrix,0); % Run lp after remove stick#1 (column=1) and strand (row=k)
        fval_sum=fval+weight_matrix(k,1) ;
               
        for j=1:strand_number-1
            p=lp_reduced(j)
            if p==0         % index before remove row k
                lp_modify(j)=0
            end
            if p~=0
                lp_modify(j)=lp_reduced(j)+1
            end
        end
        
        w=1
        for i=1:strand_number
            if i~=k
                lp_topology(i).num_strand=i;
                lp_topology(i).num_stick=lp_modify(w);
                w=w+1;
            end
            
            if i==k
                lp_topology(i).num_strand=k;
                lp_topology(i).num_stick=1;
            end
        end

        lp_reduced_matrix(k).k=k;
        lp_reduced_matrix(k).lp_topology=[lp_topology];
        lp_reduced_matrix(k).lp_score=fval_sum;       
    end
    
    % Sort topologies based on the minimum score of LP
    table_lp_matrix=struct2table(lp_reduced_matrix);
    sorted_table=sortrows(table_lp_matrix,'lp_score');
    sorted_table(:,1)=[];
    sorted_struct_lp=table2struct(sorted_table)
    
    % Determine rank for each topology
    for i=1:strand_number
        sorted_struct_lp(i).rank=i
        sorted_struct_lp(i).lp_topology=sorted_struct_lp(i).lp_topology;
        sorted_struct_lp(i).lp_score=sorted_struct_lp(i).lp_score;
    end
     
    
    all_topologies=[sorted_struct_lp.lp_score]
    [minimum,rank_topology_sheet] = min(all_topologies)
    final_topology_sheet=sorted_struct_lp(rank_topology_sheet).lp_topology
    
    
    
    % ***********    Find the direction in Topology for final assignment (DDA algorithm) *******************
    
    for i=1:strand_number
        
        if final_topology_sheet(i).num_stick ~= 0
            index=final_topology_sheet(i).num_stick;
            len=length(generated_stick_voxels(index).Mat);
            c=1;
            inverse_stick_voxels = zeros(len,3);  % preallocate of matrix
            for j=len:-1:1
                inverse_stick_voxels(c,:)=generated_stick_voxels(index).Mat(j,:)
                c=c+1;
            end
            
            A = strand_voxels(i).Mat;
            B1 = generated_stick_voxels(index).Mat;
            dtw_x1=dtw(A(:,1),B1(:,1))
            dtw_y1=dtw(A(:,2),B1(:,2))
            dtw_z1=dtw(A(:,3),B1(:,3))
            dtw_distance1 = sum(dtw_x1+dtw_y1+dtw_z1);
            
            B2 =inverse_stick_voxels
            dtw_x2=dtw(A(:,1),B2(:,1))
            dtw_y2=dtw(A(:,2),B2(:,2))
            dtw_z2=dtw(A(:,3),B2(:,3))
            dtw_distance2 = sum(dtw_x2+dtw_y2+dtw_z2);
            
            
            if dtw_distance1<=dtw_distance2 direction=1
            else direction=-1
 
            end
            
                    
        else        
            direction=0
            
        end
        final_topology_sheet(i).Direction= direction;
        
    end
    
 final_topology=final_topology_sheet;
 
end % end sheet info


%%
% *****************************    Plot3d   **************************
% ************************** plot for helices   ***********************


if helix_info==1
    figure('Name','Plot for Alpha Helices', 'Color','green');
    
    for i=1:helix_number
        
        subplot(1,2,1);
        title(' Extracted Alpha Helices from the Model');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        %     box on;
        grid on;
        plot3(helix_voxels(i).Mat(:,1),helix_voxels(i).Mat(:,2),helix_voxels(i).Mat(:,3),'linewidth',3); hold on;
        x1_middle=mean(helix_voxels(i).Mat(:,1));
        y1_middle=mean(helix_voxels(i).Mat(:,2));
        z1_middle=mean(helix_voxels(i).Mat(:,3));
        str='H%d';
        text( x1_middle,y1_middle,z1_middle,sprintf(str,i));
        
    end
    
    for i=1:stick_number_helices
        subplot(1,2,2);
        title('Extracted Alpha Sticks from the Cryo-EM Map');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        %     box on;
        grid on;
        plot3(stick_voxels_helices(i).Mat(:,1),stick_voxels_helices(i).Mat(:,2),stick_voxels_helices(i).Mat(:,3),'linewidth',12); hold on;
        x1_middle=mean(stick_voxels_helices(i).Mat(:,1));
        y1_middle=mean(stick_voxels_helices(i).Mat(:,2));
        z1_middle=mean(stick_voxels_helices(i).Mat(:,3));
        str='S%d';
        text( x1_middle,y1_middle,z1_middle,sprintf(str,i));
    end
    
end  % end of helix info for plot


% *****************************    Plot3d   **************************
% *************************** plot for sheets   ***********************
if sheet_info==1
    
    figure('Name','Plot for Beta strands','Color','yellow');
    
    for i=1:strand_number
        
        subplot(1,2,1);
        title('Extracted Beta Strands from the Model');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        %     box on;
        grid on;
        plot3(strand_voxels(i).Mat(:,1),strand_voxels(i).Mat(:,2),strand_voxels(i).Mat(:,3),'linewidth',3); hold on;
        x1_middle=mean(strand_voxels(i).Mat(:,1));
        y1_middle=mean(strand_voxels(i).Mat(:,2));
        z1_middle=mean(strand_voxels(i).Mat(:,3));
        str='ST%d';
        text( x1_middle,y1_middle,z1_middle,sprintf(str,i));
        
    end
    
    
    for i=1:stick_number_strands
        subplot(1,2,2);
        title('Extracted Beta Sticks from the Cryo-EM Map');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        %     box on;
        grid on;
        plot3(generated_stick_voxels(i).Mat(:,1),generated_stick_voxels(i).Mat(:,2),generated_stick_voxels(i).Mat(:,3),'linewidth',12); hold on;
        x1_middle=mean(generated_stick_voxels(i).Mat(:,1));
        y1_middle=mean(generated_stick_voxels(i).Mat(:,2));
        z1_middle=mean(generated_stick_voxels(i).Mat(:,3));
        str='S%d';
        text( x1_middle,y1_middle,z1_middle,sprintf(str,i));
    end
    
  
end  % end of sheet info for plot



if helix_info == 1 && sheet_info ==1

    final_topology = {final_topology_helix, final_topology_sheet};
    
end


end

