library(draw)
c = c(33,-85,34,-83.5)
drawBox(x = (c[2]+c[4])/2, y = (c[1]+c[3])/2, width = 2.5, height = 1)
plot(c(0,300), c(0,600), type="n", xlab ="", ylab ="")
r = rect(angle = 90, xleft = 100, xright = 200, ybottom = 300, ytop = 450)

point = c(120, 400)
Largeur = 20
longueur = 100
angle = 45
a = c(x = point[1]+sqrt( (-Largeur*(sin(angle))/2)^2 + (Largeur*(cos(angle))/2)^2 ), y = point[2]+sqrt( (-Largeur*(sin(angle))/2)^2 + (Largeur*(cos(angle))/2)^2 ))


plot(c(0,100), c(0,100), type="n", xlab ="", ylab ="")
P = c(40, 80)
Largeur = 10
longueur = 30
angle = 45

A = c(P[1]+0.5*Largeur*sin(angle), P[2]-0.5*Largeur*cos(angle))
B = c(P[1]-0.5*Largeur*sin(angle), P[2]+0.5*Largeur*cos(angle))
C = c(C[1]+longueur*cos(angle), C[2]-longueur*sin(angle))
D = c(B[1]+longueur*cos(angle), B[2]-longueur*sin(angle))

AB = c(-Largeur*sin(angle), Largeur*cos(angle))

lines(x = A,y = B)
lines(x = B,y = C)
lines(x = C,y = D)
lines(x = D,y = A)
