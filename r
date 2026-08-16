[1mdiff --git a/src/castingSolidificationMelting/castingSolidificationMelting.C b/src/castingSolidificationMelting/castingSolidificationMelting.C[m
[1mindex d74a4b9..fb88bae 100644[m
[1m--- a/src/castingSolidificationMelting/castingSolidificationMelting.C[m
[1m+++ b/src/castingSolidificationMelting/castingSolidificationMelting.C[m
[36m@@ -158,6 +158,16 @@[m [mvoid Foam::fv::castingSolidificationMelting::readCoeffs(const dictionary& dict)[m
     }[m
 [m
     beta_ = dict.lookup<scalar>("beta");[m
[32m+[m[32m    betaC_ = dict.lookupOrDefault<scalar>("betaC", 0.0);[m
[32m+[m
[32m+[m[32m    if (mag(betaC_) > SMALL && !compositionDependentLiquidus_)[m
[32m+[m[32m    {[m
[32m+[m[32m        FatalIOErrorInFunction(dict)[m
[32m+[m[32m            << "Solutal buoyancy requires "[m
[32m+[m[32m            << "compositionDependentLiquidus = true because the "[m
[32m+[m[32m            << "buoyancy source uses the liquid composition CarbonL."[m
[32m+[m[32m            << exit(FatalIOError);[m
[32m+[m[32m    }[m
 [m
     if (mode_ == thermoMode::lookup)[m
     {[m
[36m@@ -487,6 +497,7 @@[m [mFoam::fv::castingSolidificationMelting::castingSolidificationMelting[m
     q_(NaN),[m
     pullVelocity_(vector::zero),[m
     beta_(NaN),[m
[32m+[m[32m    betaC_(0.0),[m
     alpha1_[m
     ([m
         IOobject[m
[36m@@ -626,7 +637,19 @@[m [mvoid Foam::fv::castingSolidificationMelting::addSup[m
         const scalar alpha1c = alpha1_[celli];[m
 [m
         const scalar S = -Cu_*sqr(1.0 - alpha1c)/(pow3(alpha1c) + q_);[m
[31m-        const vector Sb = rhoRef_*g*beta_*deltaT_[i];[m
[32m+[m[32m        const scalar deltaCarbon =[m
[32m+[m[32m            compositionDependentLiquidus_[m
[32m+[m[32m        ? CarbonL_[celli] - CarbonRef_[m
[32m+[m[32m        : 0.0;[m
[32m+[m
[32m+[m[32m        const vector SbThermal =[m
[32m+[m[32m            rhoRef_*g*beta_*deltaT_[i];[m
[32m+[m
[32m+[m[32m        const vector SbSolutal =[m
[32m+[m[32m            rhoRef_*g*betaC_*deltaCarbon;[m
[32m+[m
[32m+[m[32m        const vector Sb =[m
[32m+[m[32m            SbThermal + SbSolutal;[m
 [m
         Sp[celli] += Vc*S;[m
         Su[celli] += Vc*(Sb + S*pullVelocity_);[m
[1mdiff --git a/src/castingSolidificationMelting/castingSolidificationMelting.H b/src/castingSolidificationMelting/castingSolidificationMelting.H[m
[1mindex 085d102..83591e0 100644[m
[1m--- a/src/castingSolidificationMelting/castingSolidificationMelting.H[m
[1m+++ b/src/castingSolidificationMelting/castingSolidificationMelting.H[m
[36m@@ -229,6 +229,9 @@[m [mprivate:[m
         //- Thermal expansion coefficient [1/K][m
         scalar beta_;[m
 [m
[32m+[m[32m        //- Solutal expansion coefficient [per unit carbon mass fraction][m
[32m+[m[32m        scalar betaC_;[m
[32m+[m
         //- Reference heat capacity for thermoMode::lookup[m
         scalar CpRef_;[m
 [m
