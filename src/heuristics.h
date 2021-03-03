#pragma once

const uint MIN_SEG = 50;
const uint MAX_SEG = 250;

const uint MIN_SEG1 = 70;
const uint MAX_SEG1 = 170;
const double SEG_LENGTH_PENALTY1 = 10;

const uint MIN_SEG2 = 85;
const uint MAX_SEG2 = 150;
const double SEG_LENGTH_PENALTY2 = 5;

const uint MIN_AB1 = 35;
const uint MAX_AB1 = 80;
const double AB_PENALTY1 = 5;

const uint MIN_BC1 = 8;
const uint MAX_BC1 = 45;
const double BC_PENALTY1 = 5;

const uint SHORT_QUERY = 100;

const double PERMUTED_PENALTY = 10;

const double MIN_HICONF = 20;
const double MIN_LOCONF = 10;
const double MIN_POSSIBLE = 5;

const double MIN_PSSM_SCORE = 2;
const double MIN_PSSM_PENALTY = 5;

const double MIN_C_SCORE = 2;
const double MIN_PSSM_SCORE_PERMUTED = 4;

const double REWARD_DDGGDD = 10;
const double REWARD_DDGSDD = 10;
const double REWARD_DNMSDD = 10;
const double REWARD_DNxGDN = 8;
const double PENALTY_NOT_GDD_SDD_GDN = 10;

const uint MAX_X = 10;
