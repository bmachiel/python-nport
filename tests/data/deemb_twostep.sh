#!/bin/sh

hspice deemb_twostep.sp -o /tmp
hspice deemb_twostep_open.sp -o /tmp
hspice deemb_twostep_short.sp -o /tmp
