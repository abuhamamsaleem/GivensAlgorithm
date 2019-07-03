disp("Running...")
function [Q R] = GQR(A)

global A
[n,m] = size(A) %extracting the matrix size
colus = {}  %cell array to store the columns
lowerElements = []  %array to store the lower triangular elements
indices = []  %array to store the indices of the lower triangular elements
counts = 0  %counter to index the indices array with
kLowerElements = 0 %counter to index the lower Elements array with
ct = [];
st = [];
p = [];
q = [];

for j=1 : m %loop on the columns
  for i=1 : n %loop on the rows
    colus{ j } = A(:,j);  %extracting column j of the matrix A
    colu = colus{ j };  %extracting the colu j from the cell array into a vector
    if i>j  %extracting the lower elements condition
    kLowerElements++;
    lowerElements([kLowerElements]) = colu(i);  %storing the elements
    counts++;
    indices([counts]) = i; %storing the elements indices
    p([counts]) = j  %storing the column index of the to be deleted element
    q([counts]) = i  %storing the row index of the to be deleted element
    r = sqrt((A(q(counts),p(counts))).^2 + (A(p(counts),p(counts))).^2);  %calculating the diameter
    ct([counts]) = A(p(counts),p(counts))/r  %calculating the cosine component
    st([counts]) = -(A(q(counts),p(counts))/r)  %calculating the sine component
    endif
  endfor
endfor
%%%%% ep %%%%%
    %A = [ 1 2 6 ;8 9 6 ;2 5 4]
    %n = 3 %this will be later imported from size(A)
    minus = ones(n,n);
    minus(q,p) = -1;
%%%%%%%%%%%%%%
    ep = {};  %cell for storing the unit vectors of p elements    
    [pn, pm] = size(p)  %storing the size of the p vector (elements under the diagonal)
  for t = 1 : pm
    ep{ t } = zeros(n,1);
    e = ep { t };          %%creating the unit vectors and storing them in the cell ep. 
    e(p(t),1) = 1;
    ep { t } = e
    assignin('base','ep',ep)
  endfor
%%%%%%%%%%%%%%
%%%%% eq %%%%%
    eq = {};  %cell for storing the unit vectors of q elements
    [qn, qm] = size(q)  %storing the size of the q vector (elements under the diagonal)
  for b = 1 : qm
    eq{ b } = zeros(n,1);
    e = eq { b };           %%creating the unit vectors and storing them in the cell eq.
    e(q(b),1) = 1;
    eq { b } = e
    assignin('base','eq',eq)
  endfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CALCULATING THE GIVENS MATRIX %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gc = cell(pm,1); %cell for storing the givens matrices.
    Gc { 1 }; %indexing the first givens matrix for testing purposes.
    counts_1 = 1;
  for v = 1 : pm  %looping around the givens matrices with the length of the qm or pm vector
    r = sqrt((A(q(counts_1),p(counts_1))).^2 + (A(p(counts_1),p(counts_1))).^2);  %calculating the diameter
    ct([counts_1]) = A(p(counts_1),p(counts_1))/r  %calculating the cosine component
    st([counts_1]) = -(A(q(counts_1),p(counts_1))/r)  %calculating the sine component
    Gc{ v } = eye(n) + ((ct(v)-1)*(( ep { v } * ep { v }' + eq { v } * eq { v }' ))) + (st(v)*(( ep { v } * eq { v }' - eq { v } * ep { v }' )));  
    
    A = Gc{ v }' * A;
    counts_1++;
    
    %Gc{ v };  %Givens matrix equation + indexing the matrix for testing
   assignin('base','R',A)
   assignin('base','Gc',Gc)
  endfor
  R = A
Tc = Gc { pm };
for k = numel(Gc):-1:2
   Tc = Tc * Gc{k-1};
end
  %G = G * Tc;
  G = Tc';
  Q = G;
  assignin('base','Q',Tc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endfunction

global A
%%%%%%%%%%%%%%%%%%
tic();
for mat_count = 1 : 30
  if mat_count < 10
  A = 6*rand(5,5)
  [Q R] = feval('GQR',A)
  else if 10 <= mat_count < 20
  A = randi(10,10);
  [Q R] = feval('GQR',A)
  else if mat_count >= 20
  A = randi([1 1000],20);
  [Q R] = feval('GQR',A)
endif
endif
endif
endfor
elapsed_time = toc();
assignin('base','t',elapsed_time)
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
tic();
for mat_count2 = 1 : 30
  if mat_count2 < 10
  A = 6*rand(5,5)
  [Q R] = qr(A)
  else if 10 <= mat_count2 < 20
  A = randi(10,10);
  [Q R] = qr(A)
  else if mat_count2 >= 20
  A = randi([1 1000],20);
  [Q R] = qr(A)
endif
endif
endif
endfor
elapsed_time1 = toc();
assignin('base','t1',elapsed_time1)
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%[Q R] = feval('GQR',A)
Hm=@(h) eye(length(h)) - 2/(h'*h)*(h*h');
tic();
A1 = [];
for mat_count1 = 1 : 30
  if mat_count1 < 10
  A1{mat_count1} = 6*rand(5,5)
  nh = size(A1{mat_count1})
      Qh=eye(nh);
    for kh=1:nh-1
        b=norm(A1{mat_count1}(kh:end,kh),2)*eye(nh-kh+1,1);
        h=A1{mat_count1}(kh:end,kh)-b;
        H=eye(nh);
        H(kh:end,kh:end)=Hm(h);
        A1{mat_count1}=H*A1{mat_count1};
        Qh=Qh*H;
    end;
  else if 10 <= mat_count1 < 20
  A1{mat_count1} = randi(10,10);
      nh = size(A1{mat_count1})
      Qh=eye(nh);
    for kh=1:nh-1
        b=norm(A1{mat_count1}(kh:end,kh),2)*eye(nh-kh+1,1);
        h=A1{mat_count1}(kh:end,kh)-b;
        H=eye(nh);
        H(kh:end,kh:end)=Hm(h);
        A1{mat_count1}=H*A1{mat_count1};
        Qh=Qh*H;
    end;
  else if mat_count1 >= 20
  A1{mat_count1} = randi([1 1000],20);
  A = randi([1 1000],20);
  nh = size(A1{mat_count1})
      Qh=eye(nh);
    for kh=1:nh-1
        b=norm(A1{mat_count1}(kh:end,kh),2)*eye(nh-kh+1,1);
        h=A1{mat_count1}(kh:end,kh)-b;
        H=eye(nh);
        H(kh:end,kh:end)=Hm(h);
        A1{mat_count1}=H*A1{mat_count1};
        Qh=Qh*H;
    end;
endif
endif
endif
endfor
elapsed_time2 = toc();
assignin('base','t2',elapsed_time2)
