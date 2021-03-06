#                              -*- Mode: Perl -*- 
# Verb.pm -- 
# Author          : Ulrich Pfeifer
# Created On      : Thu Feb  1 09:10:48 1996
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Sun Apr  3 11:44:28 2005
# Language        : Perl
# Update Count    : 29
# Status          : Unknown, Use with caution!

package Text::German::Verb;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(%VERB);

{
  local ($_);
  
  while (<DATA>) {
    chomp;
    ($verb, $key) = split;
    $VERB{$verb} = [split ':', $key];
  }
  close DATA;
}

sub reduce {
    my($v,$s,$e) = @_;
    my $ge = ($v.$s =~ /^ge/)?'ge':'';
    
    while (1) {		# algorithmus unklar
	#print "reduce: $s\n";
	if (defined $VERB{$s}) {
	    if ($VERB{$s}->[1]) { # 'ge' gehoert zum stamm 
		my $vg = $v;
		$vg =~ s/^ge//;
		return ($vg, $ge.$VERB{$s}->[0], $e);
	    } else {
		return ($v, $VERB{$s}->[0], $e);
	    }
	}
	last unless $e;
	$s .= substr($e,0,1);
	$e  = substr($e,1);
    }
    return undef;
}

1;
# stamm ersatz:ge gehoert zum wort
# f�hr	fahr:0
__DATA__
��	ess:0
a�	ess:0
b�ck	back:0
b�nd	bind:0
b�r	b�r:1
b�rg	berg:0
b�t	bitt:0
b�g	bieg:0
b�t	biet:0
b�k	back:0
back	back:0
band	bind:0
bar	berst:0
barg	berg:0
bat	bitt:0
bei�	bei�:0
ber	berst:0
berg	berg:0
berst	berst:0
bet	bitt:0
bi�	bei�:0
bieg	bieg:0
bier	b�r:1
biet	biet:0
bind	bind:0
bir	berst:0
birg	berg:0
biss	bei�:0
bitt	bitt:0
bl�	blas:0
blas	blas:0
bleib	bleib:0
bleich	bleich:0
blich	bleich:0
blieb	bleib:0
blies	blas:0
bog	bieg:0
bor	berst:0
borg	berg:0
bot	biet:0
br�	brat:0
br�ch	bring:0
br�t	brat:0
brach	bring:0
brann	brenn:0
brat	brat:0
brech	brech:0
brenn	brenn:0
brich	brech:0
briet	brat:0
bring	bring:0
broch	brech:0
buk	back:0
bund	bind:0
d�ch	denk:0
d�ng	ding:0
d�nk	d�nk:0
d�rb	derb:0
d�rf	d�rf:0
dach	denk:0
dang	ding:0
darb	derb:0
darf	d�rf:0
deih	deih:1
denk	denk:0
derb	derb:0
deuch	d�nk:0
dieh	deih:1
ding	ding:0
dirb	derb:0
dorb	derb:0
dr�ng	dring:0
dr�sch	dresch:0
dr�sch	dresch:0
dr�ss	drie�:0
drang	dring:0
drasch	dresch:0
dresch	dresch:0
drie�	drie�:0
dring	dring:0
drisch	dresch:0
dro�	drie�:0
drosch	dresch:0
dross	drie�:0
drung	dring:0
dung	ding:0
durf	d�rf:0
ess	ess:0
f�hl	fehl:0
f�hr	fahr:0
f�ll	fall:0
f�nd	find:0
f�ng	fang:0
f�ch	fecht:0
f�hl	fehl:0
f�hr	f�hr:0
fahl	fehl:0
fahr	fahr:0
fall	fall:0
fand	find:0
fang	fang:0
fech	fecht:0
fecht	fecht:0
fehl	fehl:0
fich	fecht:0
fiehl	fehl:0
fiel	fall:0
find	find:0
fing	fang:0
fl�ch	flecht:0
fl�g	flieg:0
fl�h	flieh:0
fl�ss	fliess:0
flech	flecht:0
flecht	flecht:0
flei�	flei�:0
fli�	flei�:0
flich	flecht:0
flicht	flecht:0
flieg	flieg:0
flieh	flieh:0
fliess	fliess:0
fliss	flei�:0
flo�	fliess:0
floch	flecht:0
flog	flieg:0
floh	flieh:0
floss	fliess:0
foch	fecht:0
fohl	fehl:0
fr��	fress:0
fr�g	frag:0
fr�r	frier:0
fra�	fress:0
frag	frag:0
fress	fress:0
fri�	fress:0
frier	frier:0
fror	frier:0
frug	frag:0
fuhr	fahr:0
fund	find:0
g��	gess:0
g�b	geb:0
g�l	gelt:0
g�nn	ginn:0
g�r	g�r:0
g�l	gelt:0
g�nn	ginn:0
g�r	g�r:0
g�ss	gie�:0
ga�	gess:0
gab	geb:0
gal	gelt:0
gang	geh:0
gann	ginn:0
geb	geb:0
geh	geh:0
gelt	gelt:0
gess	gess:0
gi�	gess:0
gib	geb:0
gie�	gie�:0
gil	gelt:0
gilt	gelt:0
ging	geh:0
ginn	ginn:0
gl�mm	glimm:0
gleich	gleich:0
gleit	gleit:0
glich	gleich:0
glimm	glimm:0
glitt	gleit:0
glomm	glimm:0
go�	gie�:0
gol	gelt:0
gonn	ginn:0
gor	g�r:0
goss	gie�:0
gr�b	grab:0
gr�b	grab:0
grab	grab:0
greif	greif:0
griff	greif:0
grub	grab:0
h�l	halt:0
h�lf	helf:0
h�lt	halt:0
h�ng	h�ng:0
h�tt	hab:0
h�b	heb:0
h�b	heb:0
h�lf	helf:0
ha	hab:0
hab	hab:0
hal	halt:0
half	helf:0
halt	halt:0
hang	h�ng:0
hat	hab:0
hau	hau:0
heb	heb:0
hei�	hei�:0
helf	helf:0
hie�	hei�:0
hieb	hau:0
hiel	halt:0
hilf	helf:0
hing	h�ng:0
hob	heb:0
holf	helf:0
hub	heb:0
i�	ess:0
ies	
k�m	komm:0
k�mm	komm:0
k�nn	k�nn:0
k�r	kies:0
kam	komm:0
kann	k�nn:0
kies	kies:0
kl�ng	kling:0
kl�mm	klimm:0
klang	kling:0
klimm	klimm:0
kling	kling:0
klomm	klimm:0
klung	kling:0
kneif	kneif:0
kniff	kneif:0
komm	komm:0
konn	k�nn:0
kor	kies:0
kr�ch	kriech:0
kreisch	kreisch:0
kriech	kriech:0
krisch	kreisch:0
kroch	kriech:0
l��	lass:0
l�d	lad:0
l�g	lieg:0
l�ng	ling:1
l�s	les:0
l�uf	lauf:0
l�g	l�g:0
l�r	lier:0
l�sch	l�sch:0
l�d	lad:0
l�g	l�g:0
lad	lad:0
lag	lieg:0
lang	ling:1
las	les:0
lass	lass:0
lauf	lauf:0
leg	lieg:0
leid	leid:0
leih	leih:0
les	les:0
lie	les:0
lie�	lass:0
lief	lauf:0
lieg	lieg:0
lieh	leih:0
lier	lier:0
lies	les
ling	ling:1
lisch	l�sch:0
litt	leid:0
log	l�g:0
lor	lier:0
losch	l�sch:0
lud	lad:0
lung	ling:1
m��	mess:0
m�ch	m�g:0
m�g	m�g:0
m�lk	melk:0
m��	m�ss:0
m�ss	m�ss:0
ma�	mess:0
mag	m�g:0
mahl	mahl:0
meid	meid:0
melk	melk:0
mess	mess:0
mi�	mess:0
mied	meid:0
milk	melk:0
moch	m�g:0
molk	melk:0
mu�	m�ss:0
n�hm	nehm:0
n�s	nes:1
n�ss	nie�:1
nahm	nehm:0
nann	nenn:0
nas	nes:1
nehm	nehm:0
nenn	nenn:0
nes	nes:1
nie�	nie�:1
nimm	nehm:0
no�	nie�:1
nomm	nehm:0
noss	nie�:1
pf�hl	pfehl:0
pf�nd	pfind:0
pf�ng	pfang:0
pf�hl	pfehl:0
pfahl	pfehl:0
pfand	pfind:0
pfang	pfang:0
pfehl	pfehl:0
pfeif	pfeif:0
pfiehl	pfehl:0
pfiff	pfeif:0
pfind	pfind:0
pfing	pfang:0
pfl�g	pfleg:0
pfleg	pfleg:0
pflog	pfleg:0
pfohl	pfehl:0
pfund	pfind:0
preis	preis:0
pries	preis:0
qu�ll	quell:0
quell	quell:0
quill	quell:0
quoll	quell:0
r�	rat:0
r�ng	ring:0
r�nn	rinn:0
r�t	rat:0
r�ch	riech:0
r�nn	rinn:0
ra	rat:0
rang	ring:0
rann	renn:0
rat	rat:0
rei�	rei�:0
reib	reib:0
reit	reit:0
renn	renn:0
ri�	rei�:0
rie	rat:0
rieb	reib:0
riech	riech:0
rief	ruf:0
riet	rat:0
ring	ring:0
rinn	rinn:0
riss	rei�:0
ritt	reit:0
roch	riech:0
ronn	rinn:0
ruf	ruf:0
rung	ring:0
s��	sitz:0
s�h	seh:0
s�ng	sing:0
s�nk	sink:0
s�nn	sinn:0
s�uf	sauf:0
s�ff	sauf:0
s�g	saug:0
s�nn	sinn:0
s�t	sied:0
s�tt	sied:0
sa�	sitz:0
sah	seh:0
salz	salz:0
sand	send:0
sang	sing:0
sank	sink:0
sann	sinn:0
sauf	sauf:0
saug	saug:0
sch�h	scheh:1
sch�b	schieb:0
sch�l	schelt:0
sch�ll	schall:0
sch�r	scher:0
sch�ss	schie�:0
sch�f	schaff:0
schaff	schaff:0
schah	scheh:1
schal	schelt:0
schall	schall:0
scheh	scheh:1
schei�	schei�:0
scheid	scheid:0
schein	schein:0
schelt	schelt:0
scher	scher:0
schi�	schei�:0
schie�	schie�:0
schieb	schieb:0
schied	scheid:0
schieh	scheh:1
schien	schein:0
schil	schelt:0
schilt	schelt:0
schind	schind:0
schiss	schei�:0
schl�f	schlaf:0
schl�g	schlag:0
schl�ng	schling:0
schl�ss	schlie�:0
schl�g	schlag:0
schlaf	schlaf:0
schlag	schlag:0
schlang	schling:0
schlei�	schlei�:0
schleich	schleich:0
schleif	schleif:0
schli�	schlei�:0
schlich	schleich:0
schlie�	schlie�:0
schlief	schlaf:0
schliff	schleif:0
schling	schling:0
schliss	schlei�:0
schlo�	schlie�:0
schloss	schlie�:0
schlug	schlag:0
schlung	schling:0
schm�lz	schmelz:0
schmei�	schmei�:0
schmelz	schmelz:0
schmi�	schmei�:0
schmilz	schmelz:0
schmiss	schmei�:0
schmolz	schmelz:0
schn�b	schnaub:0
schnaub	schnaub:0
schneid	schneid:0
schnitt	schneid:0
schnob	schnaub:0
scho�	schie�:0
schob	schieb:0
schol	schelt:0
scholl	schall:0
schor	scher:0
schoss	schie�:0
schr�k	schreck:0
schrak	schreck:0
schreck	schreck:0
schrei	schrei:0
schreib	schreib:0
schreit	schreit:0
schri	schrei:0
schrick	schreck:0
schrie	schrei:0
schrieb	schreib:0
schrit	schreit:0
schritt	schreit:0
schrock	schreck:0
schuf	schaff:0
schund	schind:0
schw�mm	schwimm:0
schw�nd	schwind:0
schw�ng	schwing:0
schw�ll	schwell:0
schw�mm	schwimm:0
schw�r	schw�r:0
schw�r	schw�r:0
schwamm	schwimm:0
schwand	schwind:0
schwang	schwing:0
schweig	schweig:0
schwell	schwell:0
schwieg	schweig:0
schwill	schwell:0
schwimm	schwimm:0
schwind	schwind:0
schwing	schwing:0
schwoll	schwell:0
schwomm	schwimm:0
schwor	schw�r:0
schwund	schwind:0
schwung	schwing:0
schwur	schw�r:0
seh	seh:0
send	send:0
sess	sitz:0
sied	sied:0
sieh	seh:0
sing	sing:0
sink	sink:0
sinn	sinn:0
sitz	sitz:0
soff	sauf:0
sog	saug:0
sonn	sinn:0
sot	sied:0
sott	sied:0
sp�nn	spinn:0
sp�nn	spinn:0
spann	spinn:0
spei	spei:0
spi	spei:0
spie	spei:0
spinn	spinn:0
splei�	splei�:0
spli�	splei�:0
spliss	splei�:0
sponn	spinn:0
spr�ch	sprech:0
spr�ng	spring:0
spr�ss	sprie�:0
sprach	sprech:0
sprang	spring:0
sprech	sprech:0
sprich	sprech:0
sprie�	sprie�:0
spring	spring:0
spro�	sprie�:0
sproch	sprech:0
spross	sprie�:0
sprung	spring:0
st�ch	stech:0
st�hl	stehl:0
st�k	steck:0
st�nd	steh:0
st�nk	stink:0
st��	sto�:0
st�b	stieb:0
st�hl	stehl:0
st�nd	steh:0
st�rb	sterb:0
stach	stech:0
stahl	stehl:0
stak	steck:0
stand	steh:0
stank	stink:0
starb	sterb:0
stech	stech:0
steck	steck:0
steh	steh:0
stehl	stehl:0
steig	steig:0
sterb	sterb:0
stich	stech:0
stie�	sto�:0
stieb	stieb:0
stieg	steig:0
stiehl	stehl:0
stink	stink:0
stirb	sterb:0
sto�	sto�:0
stob	stieb:0
stoch	stech:0
stohl	stehl:0
storb	sterb:0
streich	streich:0
streit	streit:0
strich	streich:0
strit	streit:0
stritt	streit:0
stunk	stink:0
sung	sing:0
sunk	sink:0
t�	tun:0
t�t	tun:0
ta	tun:0
tan	tun:0
tat	tun:0
tr�f	treff:0
tr�g	trag:0
tr�nk	trink:0
tr�t	tret:0
tr�ff	trief:0
tr�g	tr�g:0
tr�g	tr�g:0
traf	treff:0
trag	trag:0
trank	trink:0
trat	tret:0
treff	treff:0
treib	treib:0
tret	tret:0
trieb	treib:0
trief	trief:0
triff	treff:0
trink	trink:0
trit	tret:0
tritt	tret:0
troff	treff:0
trog	tr�g:0
trug	trag:0
trunk	trink:0
tu	tun:0
tun	tun:0
w�ch	wachs:0
w�nd	wind:0
w�nn	winn:1
w�sch	wasch:0
w�b	web:0
w�g	wieg:0
w�nn	winn:1
w��	wiss:0
w�chs	wachs:0
w�rb	werb:0
w�rd	werd:0
w�rf	werf:0
w�sch	wasch:0
wachs	wachs:0
wand	wind:0
wann	winn:1
warb	werb:0
ward	werd:0
warf	werf:0
wasch	wasch:0
web	web:0
wei�	wiss:0
weich	weich:0
weis	weis:0
werb	werb:0
werd	werd:0
werf	werf:0
wich	weich:0
wieg	wieg:0
wies	weis:0
will	woll:0
wind	wind:0
winn	winn:1
wir	werd:0
wirb	werb:0
wird	werd:0
wirf	werf:0
wiss	wiss:0
wob	web:0
wog	wieg:0
woll	woll:0
wonn	winn:1
worb	werb:0
word	werd:0
worf	werf:0
wr�ng	wring:0
wrang	wring:0
wring	wring:0
wrung	wring:0
wu�	wiss:0
wuchs	wachs:0
wund	wind:0
wurd	werd:0
wusch	wasch:0
z�g	zieh:0
zieh	zieh:0
zog	zieh:0
zw�ng	zwing:0
zwang	zwing:0
zwing	zwing:0
zwung	zwing:0
