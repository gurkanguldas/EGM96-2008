import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Class {
	double fac(double x) {
		double fac=1;
			for(int j=1 ; j<=x ; j++)
				fac*=j;
			if(fac==0)
				fac=1;
			return fac;
	}
	double fac(int x,int y) {
		double fac=1;
		if(x-y==0)
			for(int j=1 ; j<=x+y ; j++)
				fac/=Math.sqrt(j);
		else if(y==0)
			fac=1.0;
		else 
				for(int j=(x-y)+1 ; j<=x+y ; j++)
					fac/=Math.sqrt(j);
		
			return fac;
	}
	double fac(int x,int y,double b) {
		double fac=1;
		if(b==1 ) 
			fac=Math.sqrt(1.0/(x+y-1.0)/(x+y));
			else if(b==2) 
				fac=Math.sqrt((x-y-1.0)*(x-y)/(x+y-1.0)/(x+y));
			
		
			return fac;
	}
	double cos(double x) {
		int c=1;
		double a=1,fac=1;
		while(x>2*Math.PI) {
			x-=2*Math.PI;
		}
		for(int i=2; i<=100;i+=2) {
			fac=1;
			for(int j=1 ; j<=i ; j++)
				fac*=j;
			c*=-1;
			a+=c*Math.pow(x, i)/fac;
		}
		if(x==Math.PI/2)
			a=0;
		return a;}
    double sin(double x) {
		int c=-1;
		double a=0,fac=1;
		while(x>2*Math.PI) {
			x-=2*Math.PI;
		}
		for(int i=1; i<=100;i+=2) {
			fac=1;
			for(int j=1 ; j<=i ; j++)
				fac*=j;
			c*=-1;
			a+=c*Math.pow(x, i)/fac;
		}
		return a;}
	double tan(double x) {
		double a =sin(x)/cos(x);
		return a;
	}
	double cot(double x) {
		double a =cos(x)/sin(x);
		return a;
	}
	double[] BL(double saga , double yukari , double lamda0) {
		double[]BL=new double[2];
		double y=(saga-500000);
		double a=6378137.0,
			   f=1./298.257223563,
			   n=a*f/(2.*a-f*a);
		
		double b1=(2.*a-a*f)/2.*(1+Math.pow(n, 2.)/4+Math.pow(n, 4.)/64.),
			   b2=3./2.*n-27./32.*Math.pow(n, 3.)+269./512.*Math.pow(n, 5.),
			   b3=21./16.*Math.pow(n, 2.)-55/32*Math.pow(n, 4.),
			   b4=151./96.*Math.pow(n, 3)+417./128.*Math.pow(n, 5.),
			   b5=1097./512.*Math.pow(n, 4.);

		double fi0=yukari/b1+b2*sin(2.*yukari/b1)+b3*sin(4.*yukari/b1)+b4*sin(6.*yukari/b1)+b5*sin(8.*yukari/b1),
			   nu=(Math.pow(a, 2.)-Math.pow(a*(1.-f) , 2.))*Math.pow((cos(fi0)/(a*(1.-f))), 2.),
				   N=a/(1.-f)/Math.sqrt(1.+nu),t=tan(fi0);
		
		double fi=fi0
				 +t*Math.pow(y/N, 2.)*(-1-nu)/2.
				 +t*Math.pow(y/N, 4.)*(5.+3.*t*t+6.*nu-6.*t*t*nu-3.*nu*nu-9.*t*t*nu*nu)/24.
				 +t*Math.pow(y/N, 6.)*(-61.-90.*t*t-45.*Math.pow(t, 4.)-107.*nu+162.*t*t*nu+45.*Math.pow(t, 4.)*nu)/720.
				 +t*Math.pow(y/N, 8.)*(1385.+3633.*t*t+4095.*Math.pow(t, 4.)+1575.*Math.pow(t, 6)/40320.),
			   lamda=lamda0/180.*Math.PI
			   	 +y/N/cos(fi0)
			   	 +Math.pow(y/N, 3.)/6./cos(fi0)*(-1.-2.*t*t-nu)
			   	 +Math.pow(y/N, 5.)/120./cos(fi0)*(5.+28.*t*t+24.*Math.pow(t, 4.)+6.*nu+8.*t*t*nu)
			   	 +Math.pow(y/N, 7.)/5040./cos(fi0)*(-61.-662.*t*t-1320.*Math.pow(t, 4.)-720.*Math.pow(t, 6.));
		BL[0]=fi;
		BL[1]=lamda;
		return BL;
		
	}
	double[] XYZ(double B , double L , double h) {
		double[]xyz=new double[3];
		double a=6378137.0,
			   f=1./298.257223563;
		double N=a/(1-f)/Math.sqrt(((1/(1-f)/(1-f)-1)*cos(B)*cos(B)+1));
		
			   xyz[0]=(N+h)*cos(B)*cos(L);
			   xyz[1]=(N+h)*cos(B)*sin(L);
			   xyz[2]=(N*(1-f)*(1-f)+h)*sin(B);
		
		return xyz;
	}
	double[] VLR(double x , double y , double z) {
		double[] vlr = new double[3];
		vlr [0] = Math.atan (Math.sqrt(x*x+y*y)/z);
		vlr [1] = Math.atan (y/x);
		vlr [2] = Math.sqrt (x * x + y * y + z * z);
		return vlr;

	}
	
	
}
class Read {
	private double x[], y[], h[];
	Read(String uzanti) {
		File text = new File(uzanti);
		BufferedReader buffer;
		double x[], y[], h[];
		h=new double[20];
		x=new double[20];
		y=new double[20];
		int j=0;
		try {
			buffer = new BufferedReader(new FileReader(text));
			String okuma =buffer.readLine();
			
			while(okuma!=null) {
				h[j]=Double.parseDouble(okuma.substring(23,okuma.length()));
				y[j]=Double.parseDouble(okuma.substring(0,10));
				x[j]=Double.parseDouble(okuma.substring(10,23));
					j++;
				okuma= buffer.readLine();
			}
			buffer.close();
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.h=h;
		this.x=x;
		this.y=y;}
	
	public double[] geth() {
		return this.h;
	}
	public double[] getX() {
		return this.x;
	}
	public double[] getY() {
		return this.y;
	}}
class Read4D{
	ArrayList<Integer> n,m;
	ArrayList<Double> C;
	ArrayList<Double> S;
	Read4D(String uzanti){
		ArrayList<Integer> n = new ArrayList<Integer>();
		ArrayList<Integer> m = new ArrayList<Integer>();
		ArrayList<Double> C = new ArrayList<Double>();
		ArrayList<Double> S = new ArrayList<Double>();
		File text = new File(uzanti);
		try {
			BufferedReader buffer = new BufferedReader(new FileReader(text));
			String okuma =buffer.readLine();
			
			while(okuma!=null) {
				 byte k1=0,k2=0,k3=0,k4=0,k5=0;
				 byte bosluk0 = 0,bosluk1=0,bosluk2=0,bosluk3=0,bosluk4=0;
				for(int i=1 ; i<okuma.length()+1 ; i++) {
					if(okuma.substring(i-1, i).contentEquals(" ")&&k1==0) {
						bosluk0+=1;
					 }
					else if(!okuma.substring(i-1, i).contentEquals(" ")&&bosluk1==0) {
						k1+=1;
					}
					else if(okuma.substring(i-1, i).contentEquals(" ")&&k2==0) {
						bosluk1+=1;
					}
					else if(!okuma.substring(i-1, i).contentEquals(" ")&&bosluk2==0) {
						k2+=1;
					}
					else if(okuma.substring(i-1, i).contentEquals(" ")&&k3==0) {
						bosluk2+=1;
					}
					else if(!okuma.substring(i-1, i).contentEquals(" ")&&bosluk3==0) {
						k3+=1;
					}
					else if(okuma.substring(i-1, i).contentEquals(" ")&&k4==0) {
						bosluk3+=1;
					}
					else if(!okuma.substring(i-1, i).contentEquals(" ")&&bosluk4==0) {
						k4+=1;
					}
					else if(okuma.substring(i-1, i).contentEquals(" ")&&k5==0) {
						bosluk4+=1;
					}
				 }
				okuma = okuma.replace( 'D', 'e');
				n.add(Integer.valueOf(okuma.substring(bosluk0,bosluk0+k1)));
				m.add(Integer.valueOf(okuma.substring(bosluk0+k1+bosluk1,bosluk0+k1+bosluk1+k2)));
				C.add(Double.valueOf(okuma.substring(bosluk0+k1+bosluk1+k2+bosluk2,bosluk0+k1+bosluk1+k2+bosluk2+k3)));
				S.add(Double.valueOf(okuma.substring(bosluk0+k1+bosluk1+k2+bosluk2+k3+bosluk3,bosluk0+k1+bosluk1+k2+bosluk2+k3+bosluk3+k4)));
				
				okuma =buffer.readLine();
			}
			buffer.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.n=n;
		this.m=m;
		this.S=S;
		this.C=C;
	}
	ArrayList<Integer> getN(){
		return this.n;
	}
	ArrayList<Integer> getM(){
		return this.m;
	}
	ArrayList<Double> getS(){
		return this.S;
	}
	ArrayList<Double> getC(){
		return this.C;
	}
}
class Gravitation{
	private String Uzanti0="";
	private String Uzanti1="";
	public String setData(String uzanti) {
		this.Uzanti0=uzanti;
		return this.Uzanti0;
	}
	public String setEGM(String uzanti) {
		this.Uzanti1=uzanti;
		return this.Uzanti1;
	}
	public double[] Potential (int n){
		Class math = new Class();
		Read text = new Read(Uzanti0);
		Read4D egm = new Read4D(Uzanti1);
		
		ArrayList<Double> C =egm.getC();
		ArrayList<Double> S =egm.getS();
		ArrayList<Integer> N =egm.getN();
		ArrayList<Integer> M =egm.getM();
		
		double[][]VLR = new double [20][3];
		for(int i=0 ;i<text.getY().length; i++) {
		double[]cog=math.BL(text.getY()[i], text.getX()[i],33);
		double[]xyz= math.XYZ(cog[0],cog[1], text.geth()[i]);
		double[]vlr = math.VLR(xyz[0], xyz[1], xyz[2]);

		VLR[i][0]=vlr[0];
				VLR[i][1]=vlr[1];
						VLR[i][2]=vlr[2];
		}
		
		double w=7.292115e-5;
		double O[]=new double[20];
		double V[]=new double[20];
		double W[]=new double[20];
		
		for(int q=0 ; q<20 ; q++) {
			
		O[q]=1.0/2.0*w*w*Math.pow(VLR[q][2]*math.sin(VLR[q][0]), 2);	
		int size=3;
		ArrayList<Double> Pij=new ArrayList<Double>();
		ArrayList<Double> P=new ArrayList<Double>();
		Pij.add(1.0);
		P.add(1.0);
		Pij.add(math.cos(VLR[q][0]));
		P.add(Math.sqrt(3)*Pij.get(1));
		Pij.add((Math.sqrt(1-math.cos(VLR[q][0])*math.cos(VLR[q][0])))/Math.sqrt(2.0));
		P.add(Math.sqrt(6)*Pij.get(2));

		
		
		
		for(int i=2 ; i<=n ; i++) {
			for(int j=0 ; j<=i ; j++) {
				
				if(M.get(size)==0) {
					
					int ni =N.get(size);
					int mi =M.get(size);
					Pij.add((2.0*ni-1.0)/ni*math.cos(VLR[q][0])*Pij.get(size-i)-(ni-1.0)/ni*Pij.get(size-(2*i-1)));
					
					P.add(Math.sqrt(2*ni+1)*Pij.get(size));
					
					 
				}
				else if(M.get(size)!=0) {
					int ni =N.get(size);
					int mi =M.get(size);
					
					
					if(ni-2>=mi) 
					Pij.add((Pij.get(size-(2*i-1)))*(math.fac(ni, mi, 2))+  (2.0*ni-1.0)*Math.sqrt(1.0-math.cos(VLR[q][0])*math.cos(VLR[q][0]))*Pij.get(size-i-1)*math.fac(ni,mi,1));
			
					else
					Pij.add((2.0*ni-1.0)*Math.sqrt(1.0-math.cos(VLR[q][0])*math.cos(VLR[q][0]))*Pij.get(size-i-1)*math.fac(ni,mi,1));	
					
					P.add(Math.sqrt(4*ni+2)*Pij.get(size));
					
				}
				size++;
			}
		}
		
		
		size=0;	
		for(int i=0 ; i<=n ; i++) {
			for(int j=0 ; j<=i ; j++) {
				V[q]+=(C.get(size)*Math.cos(VLR[q][1]*j)+S.get(size)*Math.sin(VLR[q][1]*j))*P.get(size)*Math.pow(6378137.0/VLR[q][2], i+1);
				size++;
			}}
			
			
		V[q]=3.986004418e+14/6378137.0*V[q];
		W[q]=V[q]+O[q];
		}
		
		/** WRITE **/
	/*	BufferedWriter Wbuffer;
		try {
			Wbuffer = new BufferedWriter(new FileWriter("E:\\Desktop\\Ptext-50.txt",true));
			for(int i=0 ; i<Pij.size() ; i++) {
				//System.out.println(N.get(i)+"  "+M.get(i)+"  "+Pij.get(i)+"  "+P.get(i));
				Wbuffer.write(N.get(i)+"  "+M.get(i)+"  "+Pij.get(i)+"  "+P.get(i)+"\n");
			}
		
		
		Wbuffer.close();
			}catch (IOException e) {e.printStackTrace();}
	*/	

		
		return W;
	}
}
class ElipsoidalGravitation{
	private String Uzanti0;
	public String setData(String uzanti) {
		this.Uzanti0=uzanti;
		return this.Uzanti0;
	}
	public double[] Potential (){
	double a =6378137.0,
			   f =298.257223563,
			   w =7.292115e-5,
			   GM=3.986004418e+14;
			   
		Class math = new Class();
		Read text = new Read(Uzanti0);
		double[][]VLR = new double [20][3];
		for(int i=0 ;i<text.getY().length; i++) {
		double[]cog=math.BL(text.getY()[i], text.getX()[i],33);
		double[]xyz= math.XYZ(cog[0],cog[1], text.geth()[i]);
		double[]vlr = math.VLR(xyz[0], xyz[1], xyz[2]);

		VLR[i][0]=vlr[0];
				VLR[i][1]=vlr[1];
						VLR[i][2]=vlr[2];
		}
		double U[]= new double[20];
		double O[]=new double[20];
		for(int q=0 ; q<20 ; q++) {		
		O[q]=1.0/2.0*w*w*Math.pow(VLR[q][2]*math.sin(VLR[q][0]), 2);
		ArrayList<Double> J = new ArrayList();
		ArrayList<Double> P = new ArrayList();
		
		int n=10;
		
		double m=w*w*a*a*a*(1.0-1.0/f)/GM;
		P.add(1.0);
		P.add(math.cos(VLR[q][0]));
		for(int i=4 ; i<=n*4 ; i+=2)
		P.add((i-1.0)/i*2.0*math.cos(VLR[q][0])*P.get(i/2-1)-(i/2.0-1.0)/i*2.0*P.get(i/2-2));
		for(int i=2 ; i<=n*2 ; i+=2) {
			if(i==2) 
				J.add(2.0/3.0/f-m/3.0-1.0/3.0/f/f+2.0/21.0/f*m);
			else if(i>2) 
				J.add(Math.pow(-1.0, i/2.0+1.0)*3.0*Math.pow((2.0/f-1.0/f/f),i/2.0)/(i+1.0)/(i+3.0)*(1.0-i/2.0+5.0*i/2.0*J.get(0)/(2.0/f-1.0/f/f)));
		}	
				
		
		for(int i=1 ; i<=n ; i++) {
			U[q]+=Math.pow(a/VLR[q][2], 2*i)*J.get(i-1)*P.get(i+1);
		}
			U[q]=GM/VLR[q][2]*(1-U[q]);
			U[q]=U[q]+O[q];
	}
		return U;
		}
}
class Gravity{
	public double[] MediumGravity(String Uzanti) {
		double a =6378137.0,
			   f =1./298.257223563,
			   w =7.292115e-5,
			   GM=3.986004418e+14;
		double e  = f*(2.-f);
		double e2 = 1./(1./(2.*f-f*f)-1.);
		double m = w*w*a*a*a*(1-f)/GM;
		double Ge=GM/a/(a*(1-f))*(1-3./2.*m-3./14.*m*e2);
		double Gk=GM/a/a*(1+m+3./7.*e2*m);
		double k =a*(1.-f)*Gk/a/Ge-1;

		Class math = new Class();
		Read text = new Read(Uzanti);
		
		int length=text.getY().length;
		
		double mediumGravity[]=new double[length];
		
		for(int i=0 ;i<length; i++) {
		double[]cog=math.BL(text.getY()[i], text.getX()[i],33);
		double G0=Ge*(1.+math.sin(cog[0])*math.sin(cog[0])*k)/Math.sqrt(1.-e*math.sin(cog[0])*math.sin(cog[0]));
		mediumGravity[i] =G0*(1.-2./a*(1+f+m-2.*f*math.sin(cog[0])*math.sin(cog[0]))*text.geth()[i]+3./text.geth()[i]*text.geth()[i]/a/a);
		}
		return mediumGravity;
	}
}
