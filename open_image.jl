using ImageView, Images, ImageIO, ImageMagick, Printf, FileIO, LinearAlgebra, Arpack

imgpath = "/Users/Jens/Documents/Documents - Jens’s MacBook Pro/PhD/Num linear algebra course/Homework 1/our code/market_snapshots/market_snapshots_0001.jpg"

#basefilename="./india_driving1_frames"
# Load the image into an image object, here the first frame
#img=load(abspath(basefilename*"_0001.jpg"));
img=load(imgpath);
# Visualize it
imshow(img)
# Determine the image size
sz=size(img); szv=sz[1]*sz[2];
# Reshape the image to a vector:
R=float(red.(img));
G=float(green.(img));
B=float(blue.(img));
v=vcat(vec(R), vec(G), vec(B))



print("XX")
