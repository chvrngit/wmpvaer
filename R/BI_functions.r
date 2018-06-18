#' p_kprsBI computes the first derivative as a function of s using an analytical formula 
#'
#' @param s The value needed for the first derivative to equal log(lambdaBI)
#' @param N Number of columns in a data matrix
#' @param psqrd A vector of squared eigenvalues 
p_kprsBI<-function(s,N,psqrd) {

p1<-length(psqrd);
p2<-N-p1;
n<- N-1;

a<-Rmpfr::mpfr((n/2),120);
s<-Rmpfr::mpfr(s,120)
x<-Rmpfr::roundMpfr( (Rmpfr::mpfr(psqrd,120  )),120  );
x<-t(as.matrix(x));

#computation of the derivative of the log first term

.a1 <- a + s;
.a3 <- (.a1) - ((1 + p1)/4);
.a4 <- (p1 * .a3);
.a5 <- (.a1^(.a4));

cap_1<-p1*(( ((.a1^(.a4 - 1))*.a3) + (.a5 * log(.a1)) )/.a5)

#computation of the derivative of log(R21)

cap_m1<-Rmpfr::mpfr((matrix(0,p1,p1)),120  );

for(i in (1:p1)){
for(j in (i:p1)){

.b1 <- a + s;
.b3 <- (-.b1)^2;
.b4 <- x[1,j];
.b5 <- x[1,i];
.b7 <- .b3 + a * (1 - 4 * (s * .b4)) + s;
.b9 <- .b3 + a * (1 - 4 * (s * .b5)) + s;
.b10 <- sqrt(.b7);
.b11 <- sqrt(.b9);
.b12 <- a * .b5;
.b13 <- a * .b4;
.b14 <- 1 - 2 * (a/.b10);
.b15 <- 1 + 2 * .b1;
.b16 <- 1 - 2 * (a/.b11);
.b17 <- 1 - 2 * (.b12/.b11);
.b18 <- 1 - 2 * (.b13/.b10);
.b19 <- .b15 - 4 * .b12;
.b20 <- .b15 - 4 * .b13;
.b21 <- .b16 * .b14;
.b22 <- .b14 * .b19;
.b23 <- .b11 * .b10;

cap_m1[i,j]<-((a*(.b16 * .b20/(.b7 * .b10) + .b22/(.b9 * .b11)) - .b21/s)/s - 
    a * ((2 * (.b19/.b9) + 2 * (.b20/.b7))/.b23 + a * 
	(.b16 *(2 * (.b20 * (2 * (a/.b7) - .b14*(1/.b10	+ 2 * (.b13/(.b7 *.b18))))/.b7)
	- 4 * (.b14/(s * .b10)))/.b11 +2*(.b22 * (2 * (a/.b9) - .b16 * 
	(1/.b11 + 2*(.b12/(.b9 *.b17))))/(.b9 * .b10)))*.b5*.b4/(s*.b17*.b18)))/
	(.b21*(1-4*(a^2*.b5*.b4/(.b17*.b18*.b11*.b10)))/s + 4*(a/.b23));    

}
}

#  then sum across columns

cap_m5b<-Rmpfr::mpfr((matrix(NA,p1,1)),120  );
for(i in (1:p1)) {
cap_m5b[i,1]<-Rmpfr::mpfr((sum(cap_m1[i,(1:p1)])),120  );
}

#  and sum across rows
cap_m5c<-Rmpfr::mpfr((sum(cap_m5b[(1:p1),1])),120  );
# multiply by -1/2 and this is finished
cap_m5<-Rmpfr::mpfr((((-.5)*(cap_m5c))),120  );
#

thrdm<-Rmpfr::mpfr((matrix(NA,1,p1)),120 );

for(i in (1:p1)){


 .c1 <- a + s; 
 .c2 <- x[1,i];
 
thrd1a<- (-(a*(1+(2*.c1)-(4*a*.c2))));
thrd1b<- (2*(((-.c1)^2) + a*(1-(4*(s*.c2)) + s)));
thrd1<-thrd1a/thrd1b;


.d1 <- a + s; 
.d2 <- x[1,i]; 
.d6 <- (-.d1)^2 + a * (1 - 4 * (s * .d2)) + s; 
.d7 <- sqrt(.d6); 
#.d8 <- 1 - 2 * (a/.d7); 
.d8<-1;
#.d10 <- log(.d8) - log(s); 
.d10<-1;

thrd2<-(.d10^s) * log(.d10)
 + s*(a*(1 + (2*.d1) - 4*(a*.d2))/((.d6 *.d7) - .d8/s)*(.d10^(s - 1))/.d8);


.e1 <- x[1,i]; 
.e2 <- a + s; 
.e6 <- (-.e2)^2 + a * (1 - 4 * (s * .e1)) + s; 
.e7 <- a * .e1; 
.e8 <- sqrt(.e6);

thrd3<-( -(a^2*(1+2*.e2-4*.e7)*.e1/(.e6*(1-2*(.e7/.e8))*.e8)));


thrdm[1,i]<-thrd1+thrd2+thrd3;
}

# sum columns
cap_m6<-Rmpfr::mpfr((sum(thrdm[1,1:p1])),120  );
# finally compute capdelta_s
capdelta_s<-Rmpfr::mpfr((cap_1+cap_m5+cap_m6),120  ); 

# digamma part
dgm<-Rmpfr::mpfr((matrix(NA,1,p1)),120 );
for(i in (1:p1)){
.f1 <- (i - 1)/2; 
dgm[1,i]<-digamma(((n - p2)/2) + s- .f1) - digamma(a+s- .f1);
}
digamma_part<-sum(dgm[1,1:p1])

# put everything together
first_der<-digamma_part + capdelta_s;

return(first_der)
}
#' p_twoF1BI, computation of an approximate 2F1 value for block independence 
#'
#' @param s The value needed for the first derivative to equal log(lambdaBI) 
#' @param N number of columns in the data matrix 
#' @param psqrd The a vector of squared eigenvalues 
p_twoF1BI<- function(s,N,psqrd){

p1<-length(psqrd);
p2<-N-p1;
n<-N-1;

a<-Rmpfr::mpfr((n/2),precBits=120);
b<-a
c<-Rmpfr::mpfr((a+s),precBits=120);
mfact<-(b/(a*(c-a)))

x<-Rmpfr::roundMpfr(Rmpfr::mpfr(psqrd,precBits=120  ),120);
x<-t(as.matrix(x))


#yhat calculations

tau<-(x*(b-c))-c;
f1<-sqrt((tau^2)- ((4*a*x)*(c-b)))- tau;
yht<-(2*a)*((f1)^(-1));
yhat<-(as.matrix(yht));


#R21 calculations

m1<-matrix((1/3),p1,p1);
M1<-Rmpfr::mpfr((m1),precBits=120  )  ##-> "mpfrMatrix"

for(i in (1:p1) ){
for(j in (i:p1) ){
M1[i,j]<-(yhat[1,i] * (yhat[1,j]/a));
}
}

m2<-matrix((1/3),p1,p1);
M2<-Rmpfr::mpfr((m2),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p1)) {
for(j in (i:p1)) {
M2[i,j]<-( (1-yhat[1,i]) * ((1-yhat[1,j])/(c-a)) );
}
}

m3<-matrix((1/3),p1,p1);
M3<-Rmpfr::mpfr((m3),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p1)) {
for(j in (i:p1)) {


M3[i,j]<-((-1)*mfact*( (x[1,i]*yhat[1,i]*(1-yhat[1,j])) /
(1-(x[1,i])*yhat[1,j]))*( ( (x[1,j]*yhat[1,j]*(1-yhat[1,i])) /
(1-(x[1,j])*yhat[1,i]))) );

}
}


M4<-(M1+M2+M3);

# multiplication

m5<-matrix(NA,p1,1);
M5<-Rmpfr::mpfr((m5),precBits=120  )  ##-> "mpfrMatrix"

for (i in (1:p1)){
M5[i,1]<-abs(prod(M4[i,(1:p1)]));
}

# final answer for R21
R21<-(prod(M5[(1:p1),1]));

#final computation for 2F1 
m6a<-(yhat/a)^(a);
m6b<-((1-yhat)/(c-a))^(c-a);
m6c<-((1-(x*yhat))^(-b));
m6<-(m6a*m6b*m6c);
pm6<-prod(m6);

twoF1a<-c^((p1*c)-((p1*(p1+1)/4)));
twoF1b<-R21^(-(1/2));
twoF1<-twoF1a*twoF1b*pm6;
return(twoF1);
}

#' p_mgammaBI, computation of a multivariate gamma function needed by p_lmglmBI 
#'
#' @param p The size or length of a required eigenvector 
#' @param a A variable computed in p_lmglmBI
p_mgammaBI<- function(p,a)   {

ivec<-c(1:p);
ivecm<-Rmpfr::mpfr(ivec,120 );
ivecm1<-(1-ivecm)/2;
ivecm2<-a+ivecm1;

ivecm2<-t(as.matrix(ivecm2))

if(any(ivecm2<=0)==TRUE){
for(i in (1:ncol(ivecm2))){
if(ivecm2[1,i]<=0) {ivecm2[1,i]<-1} 
}
}

vg<-gamma(ivecm2)
mg<-prod(vg);
piv<-Rmpfr::mpfr(pi,120)
fct<-Rmpfr::mpfr((p*((p-1)/4)),120);
firstfact<-piv^fct;
g_u_pfa<- mg*firstfact;
g_u_pfa
return(g_u_pfa);
}

#' p_lmglmBI, computation of log(mgf) of the log(lambdaBI) statistic. 
#'
#' @param s_hat, The value needed for the first derivative to equal log(lambdaBI), found by function p_find_sBI
#' @param N the number of columns in the data matrix
#' @param psqrd A vector of squared eigenvalues 
p_lmglmBI<-function (s_hat,N,psqrd) {

if(is.na(s_hat)){
ans_lmglm<-NA;
} else {

p1<-length(psqrd);
p2<-N-p1;
n<-N-1;
sv1<-Rmpfr::mpfr((n/2),precBits=120 );
sv2<-Rmpfr::mpfr((((n-p2)/2)+s_hat),precBits=120 );
sv3<-Rmpfr::mpfr(((n/2)+s_hat),precBits=120  );
sv4<-Rmpfr::mpfr(((n-p2)/2),precBits=120);
factor1<-p_mgammaBI(p1,sv1);
factor2<-p_mgammaBI(p1,sv2);
factor3<-p_mgammaBI(p1,sv3);
factor4<-p_mgammaBI(p1,sv4);

factor5<-p_twoF1BI(s_hat,N,psqrd);

ans_mglm1<-factor1*factor2;
ans_mglm2<-factor3*factor4;
ans_mglm3<-ans_mglm1/ans_mglm2;
factor6<-(det(diag(,p1)-diag(psqrd,nrow=p1)))^(n/2);
ans_mglm<-abs(ans_mglm3*factor6*factor5);

ans_lmglm<-Rmpfr::mpfr((log(ans_mglm)),precBits=120  )

}
return(ans_lmglm);
}
#' p_find_sBI performs a bisection search, finds an s such that First-Der(s)=log(lambdaBI)
#'
#' @param lambdaBI The value of the Block Independence statistic 
#' @param N the number of columns in the data matrix
#' @param psqrd The a vector of squared eigenvalues  
p_find_sBI<-function(lambdaBI,N,psqrd) {

lambdaBIn<-Rmpfr::mpfr(lambdaBI,120)

difans1<-Rmpfr::mpfr((log(lambdaBIn)),120  );
difansmx1<-Rmpfr::mpfr(matrix(NA,1,2),120)

# establish s_hat_max and s_hat_min

for(s_hat in seq(from=(-10),to=(10),by=.1)){ 

difans2<-Rmpfr::mpfr(p_kprsBI(s_hat,N,psqrd),120);

if((is.nan(difans2))==TRUE){next}
if((is.infinite(difans2))==TRUE){next}
difans<-Rmpfr::mpfr((difans1-difans2),120);
difansmx2<-Rmpfr::mpfr(matrix(NA,1,2),120)
difansmx2[1,2]<-abs(difans);
difansmx2[1,1]<-s_hat;
difansmx1<-Rmpfr::rbind(difansmx1,difansmx2)
}

difansmx1<-difansmx1[-1,]
md<-min(abs(difansmx1[,2]))
reduce<-difansmx1[abs(difansmx1[,2])==md,,drop=FALSE]
s_hat_mid<-reduce[1,1];

smm<-c((.95*s_hat_mid),(1.05*s_hat_mid))
s_hat_max<-max(smm);
s_hat_min<-min(smm);


counter<-0;
difans1<-Rmpfr::mpfr((log(lambdaBI)),120  );
difans_m<-Rmpfr::mpfr(matrix(NA,1000,1),120)

repeat{

counter<-(counter+1);

difans2<-p_kprsBI(s_hat_mid,N,psqrd);
difans<-difans1-difans2;
difans_m[counter,1]<-difans;

if(counter>1){
if ((abs(difans_m[counter,1]-difans_m[counter-1,1])<.0000000000001)){
s_hat<-Rmpfr::mpfr(NA,120  );
break}
}

if ((difans<0)==TRUE ){
s_hat_max<- Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),120  );
s_hat_mid<-Rmpfr::mpfr((s_hat_min+ ((s_hat_max-s_hat_min)/2)),120  );
}

if ((difans>0)==TRUE ) {
s_hat_min<-Rmpfr::mpfr((s_hat_min+((s_hat_max-s_hat_mid)/2)),120  );
s_hat_mid<-Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),precBits=120  );
}


difmaxmin<-Rmpfr::mpfr((abs(s_hat_max-s_hat_min)),120  );
s_hat<-Rmpfr::mpfr(s_hat_mid,120  );

if (  ((abs(difans))==0)==TRUE){ break};
if (  ((abs(difans))<=.00000001)==TRUE){ break};
if (  (difmaxmin<=.000000000001)==TRUE ) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}
if (  (counter==10000)==TRUE) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}

}

ans<-Rmpfr::mpfr(s_hat,120  )

return(ans)
}

#' p_kprs2BI computes the 2nd derivative as a function of s_hat, 
#'           s_hat found by the function p_find_sBI 
#' @param s_hat The value needed for the first derivative to equal log(lambdaBI),
#'              found by p_find_s
#' @param N the number of columns in the data matrix
#' @param psqrd The a vector of squared eigenvalues  
p_kprs2BI<- function(s_hat,N,psqrd) {


if(is.na(s_hat)){
second_der<-NA;
} else {

parma<-10^(-8);
parm1<-(s_hat+parma);
parm2<-(s_hat-parma);

first_der1<-(p_kprsBI(parm1,N,psqrd));
first_der2<-(p_kprsBI(parm2,N,psqrd));
diffr<-(first_der1)-(first_der2)
second_der<-  diffr/(2*parma)
return(second_der);
}
}

#' p_l_and_rBI computes the Lugannani and Rice approx. of power 
#'
#' @param lambdaBI The value of a block independence likelihood ratio statistic 
#' @param N the number of columns in the data matrix
#' @param psqrd The a vector of squared eigenvalues  
p_l_and_rBI<- function(lambdaBI,N,psqrd) {

s_hat<-p_find_sBI(lambdaBI,N,psqrd);

if(is.na(s_hat)){
pv_answer<-NA;
} else {

lambdaBIn<-Rmpfr::mpfr(lambdaBI,precBits=120 ) 
y<-log(lambdaBIn)
k_s_hat<-p_lmglmBI(s_hat,N,psqrd);
u1<-p_kprs2BI(s_hat,N,psqrd);
u<-(s_hat)*(sqrt(abs(u1)))
r1<-y*s_hat
r2<-r1-k_s_hat;
r3<-(2*(abs(r2)))^(1/2)
r<-(sign(s_hat))*r3;
mpower1<-Rmpfr::mpfr( (Rmpfr::pnorm(r,lower.tail=FALSE)) ,precBits=120    );
piv<-Rmpfr::mpfr(pi,120)
nd_r1<-Rmpfr::mpfr((1/(sqrt(2*piv))),precBits=120)
nd_r2<-(exp(-(r^2)/2) );
mpower2<-nd_r1*nd_r2
mpower3<-( ((1/u)-(1/r)))  ; 
mpower<-mpower1+(mpower2*mpower3); 
type2error<-(1-(mpower));
if((mpower>1)==TRUE){
mpower<-1;
type2error<-0;
}
if((type2error>1)==TRUE){
mpower<-0;
type2error<-1;
}
ans<-c(type2error);
pv_answer<-as.matrix(ans);
}
return(pv_answer)

}

#' statsBIf   Computes lambdaBI statistics,p-values, and shows the block matrix definitions
#'
#' \code{statsBIf} takes a matrix of numeric data, and a vector or scalar indicating columns numbers 
#'                for a block and returns the LR statistic lambdaBI, computed p values and a matrix w. column names 
#'                and numbers 
#' @param data_matrix Is a numeric matrix
#' @param block Is a scalar or vector of column numbers corresponding to block 1
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(wmpvaer);
#'     data(plastic);
#'     plasticm<-as.matrix(plastic[,1:3]);
#'     block_v<-c(2)  
#'     statsBIf(plasticm,block_v);
#' @keywords multivariate    &   Multivariate Techniques
#' @return list with one matrix as element 
#' @export
statsBIf<-function(data_matrix,block){

data_matrix_num<-as.numeric(data_matrix);
if( (any(is.na(data_matrix_num)))==TRUE ){stop("no characters allowed in data_matrix") }

if((length(block)>floor(((ncol(data_matrix))/2)))==TRUE) { stop("#cols blk1 must be <= #cols blk2") };


N<-ncol(data_matrix);

if((is.null(colnames(data_matrix)))==TRUE){
colnames(data_matrix)<-colnames(data_matrix,do.NULL=FALSE)
};

colnumbers<-as.matrix(c(1:N));
cnamesdm<-as.matrix(colnames(data_matrix));
cnocnm<-cbind(colnumbers,cnamesdm);

if((length(block)>1)==TRUE){
blk1<-cnocnm[block,];
colnames(blk1)<-c("c#","name_blk1")
}else{
blk1<-matrix(data=(cnocnm[block,]),1,2);
colnames(blk1)<-c("c#","name_blk1");
}

if((N==2)==TRUE){
blk2<-matrix(data=(cnocnm[-block,]),1,2);
colnames(blk1)<-c("c#","name_blk1");
}else{
blk2<-cnocnm[-block,];
colnames(blk2)<-c("c#","name_blk2")
}


if((nrow(blk1)==nrow(blk2))==TRUE){
blkm<-cbind(blk1,blk2);
}

if((nrow(blk2)>nrow(blk1))==TRUE){
bd<-(nrow(blk2)-nrow(blk1));
adrow<-matrix(data=c("  ","  "), bd,2)
blk1<-rbind(blk1,adrow);
blkm<-cbind(blk1,blk2);
}

dm_cov<-cov(data_matrix);

dm1<-as.matrix(data_matrix[,block]);
dm2<-data_matrix[,-block]
sig11<-cov(dm1);
sig22<-cov(dm2);

sig12<-cov(dm1,dm2)
sig21<-cov(dm2,dm1) 

sig11_inv<-solve(sig11);
sig22_inv<-solve(sig22);

eigm<- (sig11_inv %*% sig12 %*% sig22_inv %*% sig21)
cap_p<-eigen(eigm)$values;
psqrd<-(cap_p)^2;
lambdaBI<-(det(dm_cov))/(det(sig11)*det(sig22));

pw_and_p<-p_l_and_rBI(lambdaBI,N,psqrd) 
pvalue<-Rmpfr::asNumeric((pw_and_p[1,1]));
if(is.na(pvalue)){pvalue<-c("No Answer");
lambdaBI<-format(lambdaBI,digits=8)
};
lambdaBI<-format(lambdaBI,digits=8);
pvalue<-format(pvalue,digits=8,justify=c("none"));
lanswer<-matrix(data=c(lambdaBI,pvalue),1,2);
rownames(lanswer)<-c(" ");
colnames(lanswer)<-c("lambdaBI","p value")
blkmbd<-nrow(blkm)
adrow2<-matrix(data=c("  ","  "), (blkmbd-1),2)
ansb<-rbind(lanswer,adrow2);
blkm<-cbind(blkm,ansb);

ans<-(noquote(blkm));
return(ans)

}

#' statsBIncatf Computation with block1 columns taken in combinations n at a time
#' 
#' \code{statsBIncatf} Using the function statsBIf computes stats for block independence 
#'            using columns from data_matrix at specified "n columns at a time".The program
#'            uses parallel computation if possible.  
#' @param nc_at_time A scalar that identifies the value of number of columns at a time" 
#' @param data_matrix A data matrix composed of numeric values 	
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(doParallel);
#'     library(wmpvaer);
#'     data(plastic);
#'     plasticm<-as.matrix(plastic[,1:3]);
#'     nc_at_time<-1;  
#'     statsBIncatf(plasticm,nc_at_time);       
#' @return list statsBIf matrices 
#' @export 
statsBIncatf<- function(data_matrix,nc_at_time) {
N<-ncol(data_matrix)

data_matrix_num<-as.numeric(data_matrix);
if( (any(is.na(data_matrix_num)))==TRUE ){stop("no characters allowed in data_matrix") }

if((nc_at_time>floor((N/2)))==TRUE){stop("#cols blk1 must be <= #cols blk2")};	

#section using parallel computation

'%dopar%'<-'foreach::%dopar%'

a_cores<-((parallel::detectCores( ))-1);
cm<-combn(N,nc_at_time)	;
cmn<-ncol(cm);

if(a_cores<2){
result<-list( );
for( i in 1:cmn ){
result[[i]]<-(statsBIf(data_matrix,cm[,i]));
}

} else {
EX<-c("p_find_sBI","p_kprsBI","p_kprs2BI","p_mgammaBI",
"p_lmglmBI","p_l_and_rBI","p_twoF1BI","statsBIf");
cl1<-parallel::makeCluster(a_cores);
result<-list( );
doParallel::registerDoParallel(cl1);
result<-{foreach::foreach(i=1:cmn,.export=EX) %dopar% {result[[i]]<-statsBIf(data_matrix,cm[,i])} }
parallel::stopCluster(cl1);
}
return(result)
}
