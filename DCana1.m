
%%Delete the results of previous simulation, if any
if(exist('Results.txt'))
    delete('Results.txt');
elseif(exist('Results.xls'))
    delete('Results.xls');
    if(exist('MC_values_ode.txt'))
        delete('MC_values_ode.txt');
    end
    delete('Voltages_graph.fig');
    if(exist('Currents_graph.fig'))
        delete('Currents_graph.fig');
    end
end
%%-----------------------------------------------------------------------------
prompt='Enter the file name - ';
fname=input(prompt,'s');
netlist_fileID=fopen(fname);
netlist=textscan(netlist_fileID,'%s %s %s %s %s %s %s');
fclose(netlist_fileID);
fileID1=fopen('Element_indep.txt','wt+'); %Create an empty text file Element_indep.txt for passive elements and independent sources
fileID6=fopen('Element_trs.txt','wt+'); %Create an empty text file Element_trs.txt for transistors
fileID2=fopen('VCVS.txt','wt+'); %Create an empty text file VCVS.txt for voltage controlled voltage sources
fileID3=fopen('VCCS.txt','wt+'); %Create an empty text file VCCS.txt for voltage controlled current sources
fileID4=fopen('CCVS.txt','wt+'); %Create an empty text file CCVS.txt for current controlled voltage sources
fileID5=fopen('CCCS.txt','wt+'); %Create an empty text file CCCS.txt for current controlled current sources
%%-----------------------------------------------------------------------------
%%Initialize
num_Elements=0; %Number of passive elements
num_ElementsM=0; %Number of transistors elements
num_V=0; %Number of independent voltage sources
num_I=0; %Number of independent current sources
num_Nodes=0; %Number of nodes, excluding ground (node 0)
num_VCVS=0; %Number of voltage controlled voltage sources
num_VCCS=0; %Number of voltage controlled current sources
num_CCVS=0; %Number of current controlled voltage sources
num_CCCS=0; %Number of current controlled current sources
num_Q=0; %Number of transistors
num_L=0; %Number of inductors
%%-----------------------------------------------------------------------------
for i=1:length(netlist{1})
    s=netlist{1}{i};
    switch(s(1))
        case{'R','L','C','V','I'} %For passive elements and independent sources
            fprintf(fileID1,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} '\n']);
        case{'E'} %For voltage controlled voltage sources
            fprintf(fileID2,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} ' ' ...
                netlist{6}{i} '\n']);
        case{'G'} %For voltage controlled current sources
            fprintf(fileID3,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} ' ' ...
                netlist{6}{i} '\n']);
        case{'H'} %For current controlled voltage sources
            fprintf(fileID4,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} '\n']);
        case{'F'} %For current controlled current sources
            fprintf(fileID5,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} '\n']);
        case{'M'} %For transistors
            fprintf(fileID6,[netlist{1}{i} ' ' netlist{2}{i} ' ' ...
                netlist{3}{i} ' ' netlist{4}{i} ' ' netlist{5}{i} ' ' netlist{6}{i} ' ' netlist{7}{i} '\n']);
    end
end
%%-----------------------------------------------------------------------------
%%Read the data from Element_indep.txt text file
[Name,N1,N2,value]=textread('Element_indep.txt','%s %s %s %s');
for i=1:length(Name)
    switch(Name{i}(1))
        case{'R','L','C'}
            num_Elements=num_Elements+1;
            Element(num_Elements).Name=Name{i};
            Element(num_Elements).Node1=str2num(N1{i});
            Element(num_Elements).Node2=str2num(N2{i});
            Element(num_Elements).Value=str2double(value{i});
            if(Name{i}(1)=='L')
                num_L=num_L+1;
                Inductor(num_L).Name=Name{i};
                Inductor(num_L).N1=str2num(N1{i});
                Inductor(num_L).N2=str2num(N2{i});
                Inductor(num_L).Value=str2double(value{i});
            end
        case{'V'}
            num_V=num_V+1;
            Volt_source(num_V).Name=Name{i};
            Volt_source(num_V).Node1=str2num(N1{i});
            Volt_source(num_V).Node2=str2num(N2{i});
            Volt_source(num_V).Value=str2double(value{i});
            
        case{'I'}
            num_I=num_I+1;
            Current_source(num_I).Name=Name{i};
            Current_source(num_I).Node1=str2num(N1{i});
            Current_source(num_I).Node2=str2num(N2{i});
            Current_source(num_I).Value=str2double(value{i});
    end
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end
%%-----------------------------------------------------------------------------
%%Read the data from Element_trs.txt text file
[Name,N1,N2,N3,CH,L1,W1]=textread('Element_trs.txt','%s %s %s %s %s %s %s');
for i=1:length(Name)
    switch(Name{i}(1))
        case{'M'}
            num_ElementsM=num_ElementsM+1;
            ElementM(num_ElementsM).Name=Name{i};
            des(num_ElementsM).name = Name{i};
            ElementM(num_ElementsM).Node1=str2num(N1{i});
            ElementM(num_ElementsM).Node2=str2num(N2{i});
            ElementM(num_ElementsM).Node3=str2num(N3{i});
            ElementM(num_ElementsM).Channel= CH{i};
            ElementM(num_ElementsM).L1=str2double(L1{i});
            ElementM(num_ElementsM).W1=str2double(W1{i});
    end
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),max(str2num(N3{i}),num_Nodes)));
end

%%-----------------------------------------------------------------------------

DG_num = 0;
GG_num = 0 ;
D_num = 0 ;
SE_num = 0 ;
des(1).GG = [];
des(1).DG = [];
des(1).SE = [];
des(1).D = [];
des(1).CS = [];
 l1=1;
  SE_num = 0 ;
 for g= 1:4 % getting the detiles of the circute 
    
for i=1:num_ElementsM
    DGN(i) = 0 ;
     GGN (i)= 0 ;
     betaD(i,1)=0;
    
    if (ElementM(i).Channel == 'N')
       des(i).Beta = 200;
       des(i).VTH = .7;
   else if (ElementM(i).Channel == 'P')
      des(i).Beta = 500 ; 
      des(i).VTH = 1;
       end
   end
   for j=1:num_ElementsM 
    if  (ElementM(i).Node1 ==  ElementM(j).Node3) && (i~=j) % serise 
       
       if (des(i).SE > 0) 
       des(j).SE =  des(i).SE;
       if (ElementM(i).Node3 == 0)
        ind(ElementM(i).Node1,i)= 1 ;
       else
        ind(ElementM(i).Node1,i)= 1 ;
        ind(ElementM(i).Node3,i)= -1;
       end
       else
        SE_num = SE_num + 1;
        des(i).SE = i;
        des(j).SE = i;
   if (ElementM(i).Node3 == 0)
        ind(ElementM(i).Node1,i)= 1 ;
       else
        ind(ElementM(i).Node1,i)= 1 ;
        ind(ElementM(i).Node3,i)= -1;
       end
       end
    end
     
    if (ElementM(i).Node2 ==  ElementM(j).Node2) && (ElementM(i).Node1 ~= ElementM(i).Node2)&& ...
            (i~=j) % GG connection
        if (des(i).GG > 0) 
       des(j).GG =  des(i).GG;
        else
        GG_num = GG_num + 1;
         des(i).GG = i;
         des(j).GG = i;
        end
        else if isempty(des(j).GG)
            des(j).GG = 0;
        end
    end
    if ((ElementM(i).Node1 ==  ElementM(j).Node2) || (ElementM(i).Node2 == ElementM(j).Node3))  && (i==j) ...
            && (g==1) %daiod connection
        D_num = D_num +1;
        des(i).D = 1;
    else if isempty(des(j).D)
            des(j).D = 0;
        end
    end
        
   end 
   
    % DGP = DGN + 1;
   % for j=1:num_ElementsM
    %  if (ElementM(i).Node3 ==  ElementM(j).Node2) && (ElementM(i).Node1 ~= ElementM(j).Node3) && (ElementM(i).Channel == 'P')  && (i~=j)
      
     % des(DGP,i) ="DP";
      % des(DGP,j) ="GP";
       %DGP = DGP + 1;
     % else
     %  des(DGP,i) = 0;
    %   des(DGP,j) = 0;      
   %  end
    %end  
   end

 end
        l1=1;
        l2=1;
        for i=1:num_ElementsM
            for j=1:num_ElementsM
         if (ElementM(i).Node1 ==  ElementM(j).Node2) &&(des(i).SE ~= des(j).SE ) && (des(i).GG == 0) (i~=j) %DG connection 
      if  isempty(des(j).DG) || des(j).DG == "0"
         DG_num = DG_num+1;
        des(i).DG ='D';
        des(j).DG ='G';
      else 
        des(i).DG = des(i).DG + "D" ;
        des(j).DG = des(j).DG+ "G";
      end
         else if isempty(des(j).DG)
                 des(j).DG = "0";
                 
                 
             end
         end
      if ((ElementM(i).Node1 ==  ElementM(j).Node2)|| (ElementM(i).Node3 ==  ElementM(j).Node2))&& ...
              (des(i).SE == des(j).SE ) &&(des(i).D==0) && (i~=j)&&...
             (des(j).D==0) %casscade connection 
       if (des(i).CS > 0) &&(des(j).CS > 0)
       des(j).CS =  des(i).CS;
        else
       
         des(i).CS = i;
         des(j).CS = i;
        end
        else if isempty(des(j).CS)
            des(j).CS = 0;
            end
      end
            end
        end
        
   for   i=1:num_ElementsM 
           if  (des(i).DG ~= "0")
           DGN (des(i).SE ) = DGN (des(i).SE)+ 1 ;
                else if (des(i).GG > 0) 
             GGN (des(i).SE ) = GGN (des(i).SE)+ 1 ; 
                    end  
           end
          if (betaD(des(i).SE,l1)==0)
       for j=1:num_ElementsM
       
       if (des(i).SE == des(j).SE)
           betaD(des(i).SE,l1) = des(j).Beta;
           betaD1(des(i).SE,l2) = des(j).Beta;
           VGG(i,l2) = 0 ;
           VDS(i,l2) = 0 ;
           if (des(j).D > 0)&& (ElementM(j).Channel =='N')
            DiodeDN(des(i).SE,l1) =   des(j).VTH;
            DiodeDN1(des(i).SE,l2) =   des(j).VTH;
           else
            DiodeDN(des(i).SE,l1) = 0;
            DiodeDN1(des(i).SE,l2) = 0;
           end 
           if (des(j).D > 0)&& (ElementM(j).Channel =='P')
              DiodeDP(des(i).SE,l1) =   des(j).VTH;
               DiodeDP1(des(i).SE,l2) =   des(j).VTH;  
           else 
              DiodeDP(des(i).SE,l1) =   0; 
              DiodeDP1(des(i).SE,l2) =   0;
           end
             l1=l1+1;
             l2=l2+1;
       end
      
       end
          end
       l1=1;
       end
   

% calculate totla ID and total Beta 
  [l w] = size (betaD);
  
for i=1:l
   if betaD (i,1) ~= 0
     b= 0 ;
   for t=1 : w
        if (betaD (i,t) ~= 0)
       b= (1/sqrt(betaD(i,t)))+b;
       
        end
   end
   BetaT(i) = (1/b)^2;
   ID(i) = BetaT(i)*( Volt_source(1).Value -((sum(DiodeDP (i,:)))+ (sum(DiodeDN (i,:)))))^2 ;
   end
end
%correction for incidence matrix
for   i=1:num_ElementsM 
    if (ElementM(i).Node1 ==1)
        Vddm = i;
       ind(ElementM(i).Node1,i) = -1;
       ind(ElementM(i).Node3,i) = -1;
    end
end
ind (num_ElementsM+1, num_ElementsM+1)= 0;
%{
% taking nullter cascode into consdrtion
for   i=1:num_ElementsM 
if des(i).CS ~=0
    for j=1:num_ElementsM
if (des(i).CS == des(j).CS)&& (i~=j)
    
           
    end
end
end
 
 k = find(ind(1,:) == -1);
  [lk wk] = size(k);
  [l w] = size (ind);
  inCS(w)=0;
for   i=1:wk
    j = find(ind(:,k(i)) == -1)   
      m = k(i);
    while  ~(isempty (j))
       j = find(ind(:,m) == -1);
       [lj wj] = size(j);
       if (k(i)==m)&&(des(k(i)).CS ~= 0) 
           for n=1:num_ElementsM 
           if (des(k(i)).CS == des(n).CS)&& (k(i)~=n)
           if (ElementM(k(i)).Node2 == ElementM(n).Node3)&& (ElementM(k(i)).Channel  == 'P') 
               des(k(i)).VTH=sqrt(ID(1,des(k(i)).SE)/des(n).Beta);
           else if (ElementM(k(i)).Node2 == ElementM(n).Node3)&& (ElementM(k(i)).Channel  == 'N')
               des(k(i)).VTH=(sqrt(ID(1,des(k(i)).SE)/des(n).Beta))+(sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta));
               end
           end
           if (ElementM(k(i)).Node2 == ElementM(n).Node1)&& (ElementM(k(i)).Channel  == 'P') 
            j1 = find(ind(ElementM(n).Node1,:) == -1);
         des(k(i)).VTH=sqrt(ID(1,des(k(i)).SE)/des(j1).Beta);
           else if    (ElementM(k(i)).Node2 == ElementM(n).Node1)&& (ElementM(k(i)).Channel  == 'N') 
             j1 = find(ind(ElementM(n).Node1,:) == -1);
             des(k(i)).VTH=(sqrt(ID(1,des(k(i)).SE)/des(j1).Beta))+(sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta));
               end
           end
           end
%}
           
% building VDSs and ID
Gcon = 0 ;

[l w] = size (ind);
[bl bw]=size(betaD);
 k = find(ind(1,:) == -1);
  [lk wk] = size(k);
for   i=1:wk
      j = find(ind(:,k(i)) == -1) ;  
      m = k(i);
 outp(w)=0;
  outn(w)=0;
 c=0;
     while  ~(isempty (j))
         c= c+1;
       j = find(ind(:,m) == -1);
       [lj wj] = size(j);
       % taking nullter cascode into consdrtion
       if (k(i)==m)&&(des(k(i)).CS ~= 0) 
           for n=1:num_ElementsM 
           if (des(k(i)).CS == des(n).CS)&& (k(i)~=n)
           if (ElementM(k(i)).Node2 == ElementM(n).Node3)&& (ElementM(k(i)).Channel  == 'P') 
               des(k(i)).VTH=sqrt(ID(1,des(k(i)).SE)/des(n).Beta);
           else if (ElementM(k(i)).Node2 == ElementM(n).Node3)&& (ElementM(k(i)).Channel  == 'N')
               des(k(i)).VTH=(sqrt(ID(1,des(k(i)).SE)/des(n).Beta))+(sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta));
               end
           end
           if (ElementM(k(i)).Node2 == ElementM(n).Node1)&& (ElementM(k(i)).Channel  == 'P') 
            j1 = find(ind(ElementM(n).Node1,:) == -1);
         des(k(i)).VTH=sqrt(ID(1,des(k(i)).SE)/des(j1).Beta);
           else if    (ElementM(k(i)).Node2 == ElementM(n).Node1)&& (ElementM(k(i)).Channel  == 'N') 
             j1 = find(ind(ElementM(n).Node1,:) == -1);
             des(k(i)).VTH=(sqrt(ID(1,des(k(i)).SE)/des(j1).Beta))+(sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta));
               end
           end
           end
           end
       end
     if ind (j(lj,1), m+1) ~= 0
         m1=m+1;
     if (des(m1).CS ~= 0) 
          for n=1:num_ElementsM 
           if (des(m1).CS == des(n).CS) && (m1 < n)
           if (ElementM(m1).Node2 == ElementM(n).Node3)&& (ElementM(m1).Channel  == 'P') 
               des(m1).VTH=sqrt(ID(1,des(m1).SE)/des(n).Beta);
             else if (ElementM(m1).Node2 == ElementM(n).Node3)&& (ElementM(m1).Channel  == 'N')
               des(m1).VTH=(sqrt(ID(1,des(m1).SE)/des(n).Beta))+(sqrt(ID(1,des(m1).SE)/des(m1).Beta));
                 end
           end 
               if (ElementM(m1).Node1 == ElementM(n).Node2)&& (ElementM(n).Channel  == 'P') 
                des(n).VTH=(sqrt(ID(1,des(m1).SE)/des(n).Beta))+(sqrt(ID(1,des(m1).SE)/des(m1).Beta));
                else if (ElementM(m1).Node1 == ElementM(n).Node2)&& (ElementM(m1).Channel  == 'N')
                  des(n).VTH=sqrt(ID(1,des(m1).SE)/des(m1).Beta);
                    end
               end
               
                  
            if (ElementM(m1).Node2 == ElementM(n).Node1)&& (ElementM(m1).Channel  == 'P') 
            j1 = find(ind(ElementM(n).Node1,:) == -1);
         des(m1).VTH=sqrt(ID(1,des(m1).SE)/des(j1).Beta);
          else if    (ElementM(m1).Node2 == ElementM(n).Node1)&& (ElementM(m1).Channel  == 'N') 
             j1 = find(ind(ElementM(n).Node1,:) == -1);
             des(m1).VTH=(sqrt(ID(1,des(m1).SE)/des(j1).Beta))+(sqrt(ID(1,des(m1).SE)/des(m1).Beta));
               end     
            end
           end
          end
             
     end
 end
     
           %building VDS and VTH
                if (ElementM(k(i)).Channel  == 'P') &&(k(i)==m)
                    VDS(des(k(i)).SE,k(i)) = sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta);
                    if  (des(k(i)).GG ~= 0) && (VGG(1,des(k(i)).GG)==0)
                 VGG(1,des(k(i)).GG)= Volt_source(1).Value -  VDS(des(k(i)).SE,k(i))- des(k(i)).VTH;
                    
                    else if  (des(k(i)).GG ~= 0)&& (des(k(i)).D == 0) && (VGG(1,des(k(i)).GG)~=0)
                           des(k(i)).VTH = Volt_source(1).Value - (VGG(1,des(k(i)).GG)) - VDS(des(k(i)).SE,k(i));
                            DiodeDP1(des(k(i)).SE,k(i))  =  des(k(i)).VTH;    
                        
                        else  if (des(k(i)).DG == 'D')
                            VD(des(k(i)).SE,k(i))= Volt_source(1).Value -  VDS(des(k(i)).SE,k(i));
                 else  if (des(k(i)).DG == 'G')
                         
                   VG(des(k(i).SE),k(i))= Volt_source(1).Value -  VDS(des(k(i)).SE,k(i))- des(k(i)).VTH;
                            else if  (des(k(i)).GG ~= 0)&& (des(k(i)).D == 1) && (VGG(1,des(k(i)).GG)~=0)
                        ID (1,des(k(i)).SE)= (sum(betaD(des(k(i)).SE,c:bw))-betaD(des(k(i)).SE,1))*(VGG(1,des(k(i)).GG))^2;
                        VDS(des(k(i)).SE,k(i)) = sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta);
                        des(k(i)).VTH = Volt_source(1).Value - VGG(1,des(k(i)).GG) - VDS(des(k(i)).SE,k(i));
                        DiodeDP1(des(k(i)).SE,k(i))  =  des(k(i)).VTH  ;  
                                end
                            end
                        end
                    end
                    end
                else if (ElementM(k(i)).Channel  == 'N') &&(k(i)==m)% not fully prepared
                   VDS(des(k(i)).SE,k(i)) = sqrt(ID(1,des(k(i)).SE)/des(k(i)).Beta);
                   if  (des(k(i)).GG ~= 0) && (des(k(i)).D == 1)
                 VGG(1,des(k(i)).GG)= Volt_source(1).Value; 
                    else if (des(k(i)).GG ~= 0) && (des(k(i)).D == 0)
                              VGG(des(1,k(i)).GG)= Volt_source(1).Value - des(k(i)).VTH;
                   else if (des(k(i)).DG == 'D')
                            VD(des(k(i)).SE,m)= Volt_source(1).Value;
                    end
                        end
                    end 
                    end
                    end
 if ind (j(lj,1), m+1) ~= 0
   m=m+1;
    if (ElementM(m).Channel  == 'P') 
         if  (des(m).GG ~= 0)&& (VGG(1,des(m).GG)==0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
          VGG(1,des(m).GG)= Volt_source(1).Value - des(m).VTH - sum(VDS(des(m).SE,(k(i):m)))- ...
          sum(DiodeDP1(des(m).SE,(k(i):m)))- sum(DiodeDN1(des(m).SE,(k(i):m))) ;
          else  if (des(m).DG == 'D')
                 VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                   VD = Volt_source(1).Value - sum(VDS(des(m).SE,(k(i):m)));  
                   
          else  if (des(m).DG == 'G')
                   VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);      
                   des(m).VTH = Volt_source(1).Value - sum(VDS(des(m).SE,(k(i):m)))- VD;
       else if  (des(m).GG ~= 0)&& (des(m).D == 0) && ((VGG(1,des(m).GG))~=0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(k(i):m)))- ...
          sum(DiodeDP1(des(m).SE,(k(i):m)))- sum(DiodeDN1(des(m).SE,(k(i):m))) ;
          DiodeDP1(des(m).SE,m)  =  des(m).VTH ;               
       
           else if  (des(m).GG ~= 0)&& (des(m).D == 1) && (VGG(1,des(m).GG)~=0)&& (outp(m)==0)
         ID(1,des(m).SE)= (((1/(sum(1./sqrt(betaD(des(m).SE,c:bw)))-(1/sqrt(betaD(des(m).SE,1))))))^2)*((Volt_source(1).Value-VGG(1,des(m).GG))^2);
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(k(i):m)));
           DiodeDP1(des(m).SE,m)  =  des(m).VTH ;
           outp(m) = 1;
              m=k(i);
       
               else if  (des(m).DG == "0") &&(des(m).GG == 0)&&(des(m).CS == 0) 
              VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
                   end
                   end
               end
               end 
           end
         end
    
        
    else if (ElementM(m).Channel  == 'N')
            VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
                if  (des(m).GG ~= 0)&& (VGG(1,des(m).GG)==0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
         VGG(1,des(m).GG)= Volt_source(1).Value - des(m).VTH - sum(VDS(des(m).SE,(k(i):m-1)))- ...
          sum(DiodeDP1(des(m-1).SE,(k(i):m-1)))- sum(DiodeDN1(des(m).SE,(k(i):m-1))) ;
          else  if (des(m).DG == 'D')
                 VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                   VD = Volt_source(1).Value - sum(VDS(des(m).SE,(k(i):m)))+ VDS(des(m).SE,m);  
            else  if (des(m).DG == 'G')
                       VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                  des(m).VTH = Volt_source(1).Value - sum(VDS(des(m).SE,(k(i):m)))+ VDS(des(m).SE,m) - VD ;
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(k(i):m-1)))- ...
          sum(DiodeDP1(des(m-1).SE,(k(i):m-1)))- sum(DiodeDN1(des(m).SE,(k(i):m-1))) ;
          DiodeDP1(des(m).SE,m)  =  des(m).VTH ;               
        
                    else if  (des(m).GG ~= 0)&& (des(m).D == 1) && (VGG(1,des(m).GG)~=0)&& (outn(m)==0)
         ID(1,des(m).SE)= ((1/( sum(1./sqrt(betaD1(des(m).SE,k(i):m-1)))))^2)*(Volt_source(1).Value-VGG(1,des(m).GG))^2;
      
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(k(i):m)));
           DiodeDP1(des(m).SE,m)  =  des(m).VTH ;
           outn(m) = 1;
           m=k(i);
            else if  (des(m).DG == "0")||( (des(m).D == 0)&&(des(m).GG == 0)&&(des(m).CS == 0) )
              VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
       end   
                    end    
            
                end
                 end
                end
        end   
    end
  
else if ind (j(lj,1), m-1) ~= 0
   m=m-1;
     if (ElementM(m).Channel  == 'P') 
         if  (des(m).GG ~= 0)&& (VGG(1,des(m).GG)==0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
          VGG(1,des(m).GG)= Volt_source(1).Value - des(m).VTH - sum(VDS(des(m).SE,(m:k(i))))- ...
          sum(DiodeDP1(des(m).SE,(m:k(i))))- sum(DiodeDN1(des(m).SE,(m:k(i)))) ;
        else  if (des(m).DG == 'D')
                 VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                   VD = Volt_source(1).Value - sum(VDS(des(m).SE,(m:k(i))));  
                   
          else  if (des(m).DG == 'G')
                   VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);      
                   des(m).VTH = Volt_source(1).Value - sum(VDS(des(m).SE,(m:k(i))))- VD;
       else if  (des(m).GG ~= 0)&& (des(m).D == 0) && ((VGG(1,des(m).GG))~=0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(m:k(i))))- ...
          sum(DiodeDP1(des(m).SE,(m:k(i))))- sum(DiodeDN1(des(m).SE,(m:k(i)))) ;
          DiodeDP1(des(m).SE,m)  =  des(m).VTH ;               
       
           else if  (des(m).GG ~= 0)&& (des(m).D == 1) && (VGG(1,des(m).GG)~=0)&& (outp(m)==0)
         ID(1,des(m).SE)= (sum(betaD(des(m).SE),c:bw)-betaD(des(m).SE,1))*(VGG(1,des(m).GG))^2;
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(m:k(i))));
           DiodeDP1(des(m).SE,m)  =  des(m).VTH ;
           outp(m) = 1;
              m=k(i);
       
               else if  (des(m).DG == "0")||( (des(m).D == 0)&&(des(m).GG == 0)&&(des(m).CS == 0) )
              VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
                   end
               end
               end 
           end
         end
         end
        
    else if (ElementM(m).Channel  == 'N')
                if  (des(m).GG ~= 0)&& (VGG(1,des(m).GG)==0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
         VGG(1,des(m).GG)= Volt_source(1).Value - des(m).VTH - sum(VDS(des(m).SE,(m+1:k(i))))- ...
          sum(DiodeDP1(des(m+1).SE,(m+1:k(i))))- sum(DiodeDN1(des(m).SE,(m+1:k(i)))) ;
           else  if (des(m).DG == 'D')
                 VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                   VD = Volt_source(1).Value - sum(VDS(des(m).SE,(m:k(i))))+ VDS(des(m).SE,m);  
             else  if (des(m).DG == 'G')
                       VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);   
                  des(m).VTH = Volt_source(1).Value - sum(VDS(des(m).SE,(k(i):m)))+ VDS(des(m).SE,m) - VD ;
                else if  (des(m).GG ~= 0)&& (des(m).D == 0) && ((VGG(1,des(m).GG))~=0)
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(k(i):m-1)))- ...
          sum(DiodeDP1(des(m+1).SE,(m+1:k(i))))- sum(DiodeDN1(des(m).SE,(m+1:k(i)))) ;
          DiodeDP1(des(m).SE,m)  =  des(m).VTH ;               
        
                    else if  (des(m).GG ~= 0)&& (des(m).D == 1) && (VGG(1,des(m).GG)~=0)&& (outn(m)==0)
         ID(1,des(m).SE)= (sum(betaD(des(m).SE,1:c)))*(Volt_source(1).Value-VGG(1,des(m).GG))^2;
         VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);    
       des(m).VTH = Volt_source(1).Value - VGG(1,des(m).GG) - sum(VDS(des(m).SE,(m+1:k(i))))...
           -sum(DiodeDP1(des(m+1).SE,(m+1:k(i))))- sum(DiodeDN1(des(m+1).SE,(m+1:k(i))));
           DiodeDN1(des(m).SE,m)  =  des(m).VTH ;
           outn(m) = 1;
           m=k(i);
         
            else if  (des(m).DG == "0")||( (des(m).D == 0)&&(des(m).GG == 0)&&(des(m).CS == 0) )
              VDS(des(m).SE,m) = sqrt(ID(1,des(m).SE)/des(m).Beta);
       
                end
                end   
                    end    
            
                end
                 end
                end
    
   
    end
     end
     
   
 
   
    
     end

 end
     j = find(ind(:,m) == -1);
     end
       end

%{
%calculate VDSs and VTHS
D = 0 ;

    
for   i=1:num_ElementsM 
   if des(i).DG =='D'
        D = i;
     for j=1:num_ElementsM 
       if (des(i).SE == des(j).SE ) 
     Vds(j) = sqrt(ID(des(j).SE)/ des(j).Beta);

       end  
     end
   end
end   
%}


%{

%calculate VDSs and VTHS
Dcon = 0 ;
[l w] = size (ind);
for   i=1:num_ElementsM 
    in(i) = 0;
   if des(i).DG =='D'
        Dcon = [ElementM(i).Node1 i] ;
          if (ElementM(i).Channel == 'P')
              Vd(ElementM(i).Node1,i)= Vds(i) ;
          end
   if ind(ElementM(i).Node1,i+1) ~= 0
        for j=1:(Dcon(1,1)-1)
             if ind(ElementM(i).Node1-j,Dcon(1,2)+1) ~= 0 
               Vd(ElementM(i).Node1-j,i)= (sqrt(ID(des(i).SE)/des(Dcon(1,2)+1).Beta))+ DiodeDN(des(Dcon(1,2)+1).SE,i)+ ...
                   DiodeDP(des(Dcon(1,2)+1).SE,1);
                 Dcon(1,2) = Dcon(1,2)+1;
             end
        end
   else
            for j=1:(Dcon(1,1)-1)
             if ind(ElementM(i).Node1-j,Dcon(1,2)-1) ~= 0 
              Vd(ElementM(i).Node1-j,i)= (sqrt(ID(des(i).SE)/des(Dcon(1,2)-1).Beta))+ DiodeDN(des(Dcon(1,2)+1).SE,i)+ ...
                   DiodeDP(des(Dcon(1,2)+1).SE,1);
                 Dcon(1,2) = Dcon(1,2)-1;
             end
        end   
   end
   else if (des(i).GG > 0 ) && (des(j).SE == des(i).SE)
          GGcon = [ElementM(i).Node1 i] ;
          if (ElementM(i).Channel == 'P')
              Vg(ElementM(i).Node1,i)= Vds(i)+ des(i).VTH ;
          end
           if (ind(ElementM(i).Node1,i+1) ~= 0) 

        for j=1:(GGcon(1,1)-1)
             if ind(ElementM(i).Node1-j,GGcon(1,2)+1) ~= 0 
               Vg(ElementM(i).Node1-j,i)= (sqrt(ID(des(i).SE)/des(GGcon(1,2)+1).Beta))+ DiodeDN(des(GGcon(1,2)+1).SE,i)+ ...
                   DiodeDP(des(GGcon(1,2)+1).SE,i)  ;
                 GGcon(1,2) = GGcon(1,2)+1;
             end
        end
   else
            for j=1:(GGcon(1,1)-1)
             if ind(ElementM(i).Node1-j,GGcon(1,2)-1) ~= 0 
              Vg(ElementM(i).Node1-j,i)= (sqrt(ID(des(i).SE)/des(GGcon(1,2)-1).Beta))+ DiodeDN(des(GGcon(1,2)+1).SE,i)+ ...
                   DiodeDP(des(GGcon(1,2)+1).SE,1);
                 GGcon(1,2) = GGcon(1,2)-1;
             end
        end   
   end
       end
   end
end

% adjesting VTH 

Gcon = 0 ;

[l w] = size (ind);
[ld wd] = size(Vds);
 k = find(ind(1,:) == -1)
  [lk wk] = size(k);
for   i=1:wk
   if (k(i) > wd) && (ElementM(k(i)).Channel  == 'P')
     
       if (des(k(i)).GG > 0 ) && (des(k(i)).D == 0 )
           if (ElementM(k(i)).Channel  == 'P')
           Vds(k(i)) = sqrt(ID(des(k(i)).SE)/(des(k(i)).Beta));
           des(k(i)).VTH = (sum(Vg(:,des(k(i)).GG))) - Vds(k(i));
           else
              Vds(k(i)) = sqrt(ID(des(k(i)).SE)/(des(k(i)).Beta));
           des(k(i)).VTH = (sum(Vg(:,des(k(i)).GG))) ;   
           end
       else if (des(k(i)).DG == 'G' ) && (des(k(i)).D == 0 )
                if (ElementM(k(i)).Channel  == 'P')
       Vds(k(i)) = sqrt(ID(des(k(i)).SE))/(des(k(i)).Beta);
           des(k(i)).VTH = (sum(Vd(:,2))) - Vds(k(i));
                else 
                   Vds(k(i)) = sqrt(ID(des(k(i)).SE))/(des(k(i)).Beta);
           des(k(i)).VTH = (sum(Vd(:,2))) ;     
                end
           else if (des(k(i)).D ~= 0 )&&(des(k(i)).GG > 0 )
                
             ID(des(k(i)).SE) =    (des(k(i)).Beta) * (sum(Vg(:,des(k(i)).GG)))^2 ;
             Vds(k(i)) = sqrt(ID(des(k(i)).SE))/(des(k(i)).Beta);
         
             else if (des(k(i)).D ~= 0 )&&(des(k(i)).GG == 0 ) && (des(k(i)).DG == "0" )
                  if (ElementM(k(i)).Channel  == 'P')
                     Vds(k(i)) = sqrt(ID(des(k(i)).SE))/(des(k(i)).Beta);
                  else
                       Vds(k(i)) = sqrt(ID(des(k(i)).SE))/(des(k(i)).Beta);  
                  end
               end
           end
       end
  
   end
     j = find(ind(:,k(i)) == -1)   
     m = k(i);
     while  ~(isempty (j)) % need to make this works for both P and N 
       
       j = find(ind(:,m) == -1);
       [lj wj] = size(j);
 if ind (j(lj,1), m+1) ~= 0
     m= m +1 ;
       if (des(m).GG > 0 ) && (des(m).D == 0 )
           Vds(m) = sqrt(ID(des(m).SE)/(des(k(i)).Beta));
           Vds1(m) = sqrt(ID(des(m).SE))*sum(1./sqrt(betaD1(des(k(i)).SE,k(i):m)));
           des(m).VTH = (sum(Vg(:,des(m).GG))) - sum(Vds1(k(i):m))
       else if (des(m).DG == 'G' ) && (des(m).D == 0 )
        Vds(m) = sqrt(ID(des(m).SE)/(des(k(i)).Beta));
            Vds1(m) = sqrt(ID(des(m).SE))*sum(1./sqrt(betaD1(des(k(i)).SE,k(i):m)));
           des(m).VTH = (sum(Vd(:,2))) - sum(Vds(k(i):m));
           else if (des(m).DG == "0" ) && (des(m).GG == 0 )
               Vds(m) = sqrt(ID(des(m).SE)/(des(k(i)).Beta));      
           else if  (des(k(i)).D ~= 0 )&&(des(k(i)).GG > 0 )&&( in(m) == 0)
               in(m) = 1;   
             ID(des(m).SE) =   ( sum(1./sqrt(betaD1(des(k(i))).SE,k(i):m))) * sum(Vg(:,des(m).GG))^2 ;  
             m = k(i);
            
               end
           end
       end
 end
 
 
 end
 
     j = find(ind(:,m) == -1);
     end
   end
end

% taking R into considration 
for i=1:num_Elements
 if (Element(i).Node2==0) && (Element(i).Node1 > 0)   
    for j=1:num_ElementsM
       if ElementM(j).Node1 == Element(i).Node1
       if  des(j).D >0
           IR(i)=((Vds(j) + des(j).VTH ) /  Element(i).Value)*1000000;
           Vds(j)= sqrt(ID(j)- IR(i)/des(j).Beta);
       else 
       IR(i)=((Vds(j)  ) /  Element(i).Value)*1000000;
           Vds(j)= sqrt((ID(j)- IR(i))/des(j).Beta);
end
       
       end
    end
 end
end
%}