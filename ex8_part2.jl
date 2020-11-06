using ImageView, Images, ImageIO, ImageMagick, Printf, FileIO, LinearAlgebra, Arpack, SparseArrays
include("lanczos.jl")

basefilename="./india_driving1_frames"
basefilename = "/Users/Jens/Documents/Documents - Jens’s MacBook Pro/PhD/Num linear algebra course/Homework 1/our code/market_snapshots/market_snapshots_"

sz = [684, 1258]
szv = 860472
#m=295  # Number of frames to load
m=295
A=zeros(szv*3,m) # Matrix to store all the frames in
for k=1:m
    fname=@sprintf("%s%04d.jpg",basefilename,k);
    img=load(fname);
    #sz=size(img);
    print(size(img))
    #szv=sz[1]*sz[2];
    println(fname)
    R=float(red.(img));
    G=float(green.(img));
    B=float(blue.(img));
    va=vcat(vec(R), vec(G), vec(B));
    A[:,k]=va;
end


n = size(A)[1]
m = size(A)[2]
nm = n+m
C = spzeros(nm,nm)
C[1:n,n+1:nm] = A
C[n+1:nm,1:n] = A'

m = 10
b=randn(size(C)[1]);
b = b/norm(b)


H, Q = lanczos(C, b, m)

F = eigen(H)
svec = Q[:,1:end-1]*F.vectors # Ritz vectors
svecLast = svec[:,end]

sv = F.values[end]
uu = svecLast[1:n]
uu = uu/norm(uu)
vv = svecLast[n+1:end]
vv = vv/norm(vv)

lanApprox = uu*sv*vv'

# Take first column (representing first image)
vx = lanApprox[:,1]

vv=reshape(vx,szv,3);
R=reshape(vv[:,1],sz[1],sz[2]);
G=reshape(vv[:,2],sz[1],sz[2]);
B=reshape(vv[:,3],sz[1],sz[2]);
newimg=RGB.(R,G,B);
#imshow(newimg);
