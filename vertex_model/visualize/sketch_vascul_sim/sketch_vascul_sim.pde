String[] lines, lines2, lines3, lines4,lines5, pieces, pieces2, pieces3, pieces4,pieces5;
int i,j,jmax,k;
int ifac = 1;
int cell_type=0;
int cell_lineage=0;
int cell_defect=0;

//float a=8,b=220,c=210,d=2;//initial3,4,5,8,9,10
float a=8,b=220,c=240,d=2;//initial1,2,6

float g=10000000;
float af=1000000, bf=100;
FloatList x = new FloatList ();
FloatList y = new FloatList ();
FloatList Fx = new FloatList ();
FloatList Fy = new FloatList ();
float x_ch=0;
float y_ch=0;
float S1=0;
float S2=0;
float Theta1=0;
float Theta2=0;

float Stress_ave=0;
float Stress_vari=0;
float Anisotoropy_ave=0;
float Va_cellnumber=0;

void setup() {
  size(850,850);
  background(255);
  strokeWeight(3);
  frameRate(100);
  smooth();
}

void draw() {
  if(ifac<2532/*2782*//*2532*//*3282*/) {
   background(255);
   strokeWeight(3);
 
   lines  = loadStrings("filename/z_cellform"+ifac+".dat");
   lines3 = loadStrings("filename/z_celltype"+ifac+".dat");
   lines4 = loadStrings("filename/z_stresstensor"+ifac+".dat");
   lines5 = loadStrings("filename/z_centroid"+ifac+".dat");
 
  for(i=0; i<lines.length; i++) 
  {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++)
    {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
    
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
   cell_lineage=int(pieces3[1]);
   cell_defect=int(pieces3[2]);
   
    pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
     Theta1=float(pieces4[2]);
     Theta2=float(pieces4[3]);
     
    stroke(0);
    if (cell_type==1) {
    fill(255,255,0);//yellow-xylem
     }
    if (cell_type==2) {
    fill(12,0,200);//blue-pse
     }
   if (cell_type==7) {
    fill(12,0,200);//blue-pse-ln
     }
    if (cell_type==3) {
    fill(12,150,220);//aqua-ipc
     }
    if (cell_type==4) {
    fill(128,128,128);//gray-peri
     }
    if (cell_type==5) {
    fill(202,204,178);//lightgold-end
     }
    if (cell_type==6) {
   fill(150,12,220);//purple-opc
     }
     
    //OPC defect
   if(cell_defect==1){
      fill(100,0,0);
    }
   
 strokeWeight(3);
  beginShape();
  for (j=0; j<x.size(); j++)
  {
  vertex((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++)
  {
  x.remove(0);
  y.remove(0);
  }
  
  //centroid
     pieces5 = split(lines5[i],  ' ');
     x_ch=(b+(float(pieces5[0]))*a);
     y_ch=(c+(float(pieces5[1]))*a);
     
     //stress anisotropy
      stroke(220);
      strokeWeight(4);
      if (cell_type==1 || cell_type==2 || cell_type==3 || cell_type==7 || cell_type==6) {
      if (S1-S2>Anisotoropy_ave*0.3) {
       //line(x_ch-(S1-S2)*cos(Theta1)*0.5*g,y_ch-(S1-S2)*sin(Theta1)*0.5*g,x_ch+(S1-S2)*cos(Theta1)*0.5*30,y_ch+(S1-S2)*sin(Theta1)*0.5*g);  
      }
      }
      
     if (ifac==2781/*2781*//*2531*//*3281*/) {
        if (cell_type==1) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
        }
      }
      if (ifac==2781/*2781*//*2531*//*3281*/) {
        if (cell_type==3 || cell_type==6) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
        }
      }
      if (ifac==2531/*2781*//*2531*//*3281*/) {
        if (cell_type==2) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
        }
      }
      
   fill(0);
  textSize(15);
     //number_character(cell_type,x_ch-10,y_ch);
   //number_character(i,x_ch-10,y_ch);
  } //i_cell
  
  //pse_edge line red
   for(i=0; i<lines.length; i++) {
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
    if (cell_type==2) {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++) {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
   stroke(220,0,0);
   noFill();
 
  beginShape();
  for (j=0; j<x.size(); j++) {
  vertex((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++) {
    x.remove(0);
    y.remove(0);
    }
    }
   }//pse_edge line red
    ifac+=10;
   //saveFrame("foldername/z_form_####.png");
  println(ifac);
  } //ifac

   if (ifac==61 /*2781*//*2531*//*3281*/) {
    //saveFrame("foldername/z_form_1.png");
     }

} //void draw

void number_character(int i,float x, float y){
 text(i, x,y);
}

void mousePressed() {
  loop();  // Holding down the mouse activates looping
}

void mouseReleased() {
  noLoop();  // Releasing the mouse stops looping draw()
}
