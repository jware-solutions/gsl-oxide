//! GSL 2.8 CDF reference cases with original tolerances.
// Ported from GSL 2.8 cdf/test.c; values and tolerances retained.
// Copyright (C) 2002 Jason H. Stover and GSL contributors. GPL-3.0-or-later.
// The upstream #if 0 case for bug #39057 remains excluded, as in GSL's test run.
#![allow(clippy::excessive_precision, clippy::approx_constant)]
mod support;
use gsl_oxide_cdf::*;
use support::*;

/// Checks ugaussian.
#[test]
fn ugaussian() {
    check!(standard_normal_cdf, (0.0), 0.5, TEST_TOL0);
    check!(standard_normal_cdf, (1e-32), 0.5, TEST_TOL0);
    check!(
        standard_normal_cdf,
        (1e-16),
        0.5000000000000000398942280401,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (1e-8),
        0.5000000039894228040143267129,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (0.5),
        0.6914624612740131036377046105,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (0.7),
        0.7580363477769269852506495717,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (5.0),
        0.9999997133484281208060883262,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (10.0),
        0.9999999999999999999999923801,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (30.0),
        1.000000000000000000000000000,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (40.0),
        1.000000000000000000000000000,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (1e10),
        1.000000000000000000000000000,
        TEST_TOL0
    );
    check!(standard_normal_cdf, (-1e-32), 0.5, TEST_TOL0);
    check!(
        standard_normal_cdf,
        (-1e-16),
        0.4999999999999999601057719598,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (-1e-8),
        0.4999999960105771959856732870,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (-0.5),
        0.3085375387259868963622953894,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (-0.7),
        0.2419636522230730147493504282,
        TEST_TOL0
    );
    check!(
        standard_normal_cdf,
        (-5.0),
        0.0000002866515718791939116737523329,
        TEST_TOL1
    );
    check!(
        standard_normal_cdf,
        (-10.0),
        7.619853024160526065973343257e-24,
        TEST_TOL3
    );
    check!(
        standard_normal_cdf,
        (-30.0),
        4.906713927148187059533809288e-198,
        TEST_TOL3
    );
    check!(standard_normal_cdf, (-1e10), 0.0, 0.0);
    check!(standard_normal_survival, (0.0), 0.5, TEST_TOL0);
    check!(standard_normal_survival, (1e-32), 0.5, TEST_TOL0);
    check!(
        standard_normal_survival,
        (1e-16),
        0.4999999999999999601057719598,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (1e-8),
        0.4999999960105771959856732870,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (0.5),
        0.3085375387259868963622953894,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (0.7),
        0.2419636522230730147493504282,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (5.0),
        0.0000002866515718791939116737523329,
        TEST_TOL3
    );
    check!(
        standard_normal_survival,
        (10.0),
        7.619853024160526065973343257e-24,
        TEST_TOL3
    );
    check!(
        standard_normal_survival,
        (30.0),
        4.906713927148187059533809288e-198,
        TEST_TOL3
    );
    check!(standard_normal_survival, (1e10), 0.0, 0.0);
    check!(standard_normal_survival, (-1e-32), 0.5, TEST_TOL0);
    check!(
        standard_normal_survival,
        (-1e-16),
        0.5000000000000000398942280401,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-1e-8),
        0.5000000039894228040143267129,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-0.5),
        0.6914624612740131036377046105,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-0.7),
        0.7580363477769269852506495717,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-5.0),
        0.9999997133484281208060883262,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-10.0),
        0.9999999999999999999999923801,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-30.0),
        1.000000000000000000000000000,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-40.0),
        1.000000000000000000000000000,
        TEST_TOL0
    );
    check!(
        standard_normal_survival,
        (-1e10),
        1.000000000000000000000000000,
        TEST_TOL0
    );
}

/// Checks ugaussianinv.
#[test]
fn ugaussianinv() {
    check!(standard_normal_quantile, (0.9999997133), 5.0, 1e-4);
    check!(standard_normal_quantile, (0.9999683288), 4.0, 1e-6);
    check!(standard_normal_quantile, (0.9986501020), 3.0, 1e-8);
    check!(standard_normal_quantile, (0.977249868051821), 2.0, 1e-14);
    check!(
        standard_normal_quantile,
        (0.841344746068543),
        1.0,
        TEST_TOL3
    );
    check!(
        standard_normal_quantile,
        (0.691462461274013),
        0.5,
        TEST_TOL2
    );
    check!(
        standard_normal_quantile,
        (0.655421741610324),
        0.4,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (0.617911422188953),
        0.3,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (0.579259709439103),
        0.2,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (0.539827837277029),
        0.1,
        TEST_TOL1
    );
    check!(standard_normal_quantile, (0.5), 0.0, TEST_TOL0);
    check!(
        standard_normal_quantile,
        (4.60172162722971e-1),
        -0.1,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (4.20740290560897e-1),
        -0.2,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (3.82088577811047e-1),
        -0.3,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (3.44578258389676e-1),
        -0.4,
        TEST_TOL1
    );
    check!(
        standard_normal_quantile,
        (3.08537538725987e-1),
        -0.5,
        TEST_TOL2
    );
    check!(
        standard_normal_quantile,
        (1.58655253931457e-1),
        -1.0,
        TEST_TOL3
    );
    check!(standard_normal_quantile, (2.2750131948179e-2), -2.0, 1e-14);
    check!(standard_normal_quantile, (1.349898e-3), -3.0, 1e-8);
    check!(standard_normal_quantile, (3.16712e-5), -4.0, 1e-6);
    check!(standard_normal_quantile, (2.86648e-7), -5.0, 1e-4);
    check!(
        standard_normal_quantile,
        (7.61985302416052e-24),
        -10.0,
        1e-4
    );
    check!(
        standard_normal_inverse_survival,
        (7.61985302416052e-24),
        10.0,
        1e-4
    );
    check!(standard_normal_inverse_survival, (2.86648e-7), 5.0, 1e-4);
    check!(standard_normal_inverse_survival, (3.16712e-5), 4.0, 1e-6);
    check!(standard_normal_inverse_survival, (1.349898e-3), 3.0, 1e-8);
    check!(
        standard_normal_inverse_survival,
        (2.2750131948179e-2),
        2.0,
        1e-14
    );
    check!(
        standard_normal_inverse_survival,
        (1.58655253931457e-1),
        1.0,
        TEST_TOL3
    );
    check!(
        standard_normal_inverse_survival,
        (3.08537538725987e-1),
        0.5,
        TEST_TOL2
    );
    check!(
        standard_normal_inverse_survival,
        (3.44578258389676e-1),
        0.4,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (3.82088577811047e-1),
        0.3,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (4.20740290560897e-1),
        0.2,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (4.60172162722971e-1),
        0.1,
        TEST_TOL1
    );
    check!(standard_normal_inverse_survival, (0.5), 0.0, TEST_TOL0);
    check!(
        standard_normal_inverse_survival,
        (0.539827837277029),
        -0.1,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (0.579259709439103),
        -0.2,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (0.617911422188953),
        -0.3,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (0.655421741610324),
        -0.4,
        TEST_TOL1
    );
    check!(
        standard_normal_inverse_survival,
        (0.691462461274013),
        -0.5,
        TEST_TOL2
    );
    check!(
        standard_normal_inverse_survival,
        (0.841344746068543),
        -1.0,
        TEST_TOL3
    );
    check!(
        standard_normal_inverse_survival,
        (0.977249868051821),
        -2.0,
        1e-14
    );
    check!(standard_normal_inverse_survival, (0.9986501020), -3.0, 1e-8);
    check!(standard_normal_inverse_survival, (0.9999683288), -4.0, 1e-6);
    check!(standard_normal_inverse_survival, (0.9999997133), -5.0, 1e-4);
}

/// Checks exponential.
#[test]
fn exponential() {
    check!(
        exponential_cdf,
        (0.1, 0.7),
        1.33122100249818372e-1,
        TEST_TOL0
    );
    check!(
        exponential_cdf,
        (1e-32, 0.7),
        1.42857142857142857e-32,
        TEST_TOL0
    );
    check!(exponential_cdf, (1000.0, 0.7), 1.0, TEST_TOL6);
    check!(
        exponential_survival,
        (0.1, 0.7),
        8.66877899750181628e-1,
        TEST_TOL0
    );
    check!(exponential_survival, (1e-32, 0.7), 1.0, TEST_TOL0);
    check!(exponential_survival, (1000.0, 0.7), 0.0, TEST_TOL6);
}

/// Checks exponentialinv.
#[test]
fn exponentialinv() {
    check!(
        exponential_quantile,
        (0.13, 0.7),
        9.74834471334553546e-2,
        TEST_TOL0
    );
    check!(
        exponential_quantile,
        (1.42e-32, 0.7),
        9.94000000000000000e-33,
        TEST_TOL0
    );
    check!(
        exponential_inverse_survival,
        (0.86, 0.7),
        1.05576022814208545e-1,
        TEST_TOL0
    );
    check!(
        exponential_inverse_survival,
        (0.99999, 0.7),
        7.00003500023333508e-6,
        TEST_TOL6
    );
}

/// Checks exppow.
#[test]
fn exppow() {
    check!(exponential_power_cdf, (-1000.0, 0.7, 1.8), 0.0, TEST_TOL6);
    check!(
        exponential_power_cdf,
        (-0.1, 0.7, 1.8),
        0.4205349082867515493458053850,
        TEST_TOL0
    );
    check!(
        exponential_power_cdf,
        (-1e-32, 0.7, 1.8),
        0.4999999999999999999999999999,
        TEST_TOL0
    );
    check!(
        exponential_power_cdf,
        (0.1, 0.7, 1.8),
        0.5794650917132484506541946149,
        TEST_TOL0
    );
    check!(exponential_power_cdf, (1e-32, 0.7, 1.8), 0.5, TEST_TOL0);
    check!(
        exponential_power_cdf,
        (1000.0, 0.7, 1.8),
        0.9999999999999999999999956212,
        TEST_TOL6
    );
    check!(
        exponential_power_survival,
        (-1000.0, 0.7, 1.8),
        0.9999999999999999999999956212,
        TEST_TOL6
    );
    check!(
        exponential_power_survival,
        (-0.1, 0.7, 1.8),
        0.5794650917132484506541946149,
        TEST_TOL0
    );
    check!(
        exponential_power_survival,
        (-1e-32, 0.7, 1.8),
        0.5,
        TEST_TOL0
    );
    check!(
        exponential_power_survival,
        (0.1, 0.7, 1.8),
        0.4205349082867515493458053850,
        TEST_TOL0
    );
    check!(
        exponential_power_survival,
        (1e-32, 0.7, 1.8),
        0.4999999999999999999999999999,
        TEST_TOL0
    );
    check!(
        exponential_power_survival,
        (1000.0, 0.7, 1.8),
        0.0,
        TEST_TOL6
    );
}

/// Checks tdist.
#[test]
fn tdist() {
    check!(student_t_cdf, (0.0, 1.0), 0.5, TEST_TOL6);
    check!(student_t_cdf, (1e-100, 1.0), 0.5, TEST_TOL6);
    check!(
        student_t_cdf,
        (0.001, 1.0),
        5.00318309780080559e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.01, 1.0),
        5.03182992764908255e-1,
        TEST_TOL6
    );
    check!(student_t_cdf, (0.1, 1.0), 5.31725517430553569e-1, TEST_TOL6);
    check!(
        student_t_cdf,
        (0.325, 1.0),
        6.00023120032852123e-1,
        TEST_TOL6
    );
    check!(student_t_cdf, (1.0, 1.0), 0.75000000000000000e0, TEST_TOL6);
    check!(student_t_cdf, (1.5, 1.0), 8.12832958189001183e-1, TEST_TOL6);
    check!(student_t_cdf, (2.0, 1.0), 8.52416382349566726e-1, TEST_TOL6);
    check!(
        student_t_cdf,
        (10.0, 1.0),
        9.68274482569446430e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (20.0, 1.0),
        9.84097748743823625e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (100.0, 1.0),
        9.96817007235091745e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1000.0, 1.0),
        9.99681690219919441e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (10000.0, 1.0),
        9.99968169011487724e-1,
        TEST_TOL6
    );
    check!(student_t_survival, (0.0, 1.0), 0.5, TEST_TOL6);
    check!(student_t_survival, (1e-100, 1.0), 0.5, TEST_TOL6);
    check!(
        student_t_survival,
        (0.001, 1.0),
        4.99681690219919441e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.01, 1.0),
        4.96817007235091745e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.1, 1.0),
        4.68274482569446430e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.325, 1.0),
        3.99976879967147876e-1,
        TEST_TOL6
    );
    check!(student_t_survival, (1.0, 1.0), 2.5e-1, TEST_TOL6);
    check!(
        student_t_survival,
        (1.5, 1.0),
        1.87167041810998816e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (2.0, 1.0),
        1.47583617650433274e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (10.0, 1.0),
        3.17255174305535695e-2,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (20.0, 1.0),
        1.59022512561763752e-2,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (100.0, 1.0),
        3.18299276490825515e-3,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1000.0, 1.0),
        3.18309780080558939e-4,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (10000.0, 1.0),
        3.18309885122757724e-5,
        TEST_TOL6
    );
    check!(student_t_cdf, (-1e-100, 1.0), 0.5, TEST_TOL6);
    check!(
        student_t_cdf,
        (-0.001, 1.0),
        4.99681690219919441e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.01, 1.0),
        4.96817007235091744e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.1, 1.0),
        4.68274482569446430e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.325, 1.0),
        3.99976879967147876e-1,
        TEST_TOL6
    );
    check!(student_t_cdf, (-1.0, 1.0), 0.25, TEST_TOL6);
    check!(
        student_t_cdf,
        (-1.5, 1.0),
        1.87167041810998816e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-2.0, 1.0),
        1.47583617650433274e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-10.0, 1.0),
        3.17255174305535695e-2,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-20.0, 1.0),
        1.59022512561763751e-2,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-100.0, 1.0),
        3.18299276490825514e-3,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1000.0, 1.0),
        3.18309780080558938e-4,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-10000.0, 1.0),
        3.18309885122757724e-5,
        TEST_TOL6
    );
    check!(student_t_survival, (-1e-100, 1.0), 0.5, TEST_TOL6);
    check!(
        student_t_survival,
        (-0.001, 1.0),
        5.00318309780080559e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.01, 1.0),
        5.03182992764908255e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.1, 1.0),
        5.31725517430553570e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.325, 1.0),
        6.00023120032852124e-1,
        TEST_TOL6
    );
    check!(student_t_survival, (-1.0, 1.0), 7.5e-1, TEST_TOL6);
    check!(
        student_t_survival,
        (-1.5, 1.0),
        8.12832958189001184e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-2.0, 1.0),
        8.52416382349566726e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-10.0, 1.0),
        9.68274482569446430e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-20.0, 1.0),
        9.84097748743823625e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-100.0, 1.0),
        9.96817007235091745e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1000.0, 1.0),
        9.99681690219919441e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-10000.0, 1.0),
        9.99968169011487724e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.0, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1e-100, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.001, 2.0),
        5.00353553302204959e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.01, 2.0),
        5.03535445520899514e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.1, 2.0),
        5.35267280792929913e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.325, 2.0),
        6.11985772746873767e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1.0, 2.0),
        7.88675134594812882e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1.5, 2.0),
        8.63803437554499460e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (2.0, 2.0),
        9.08248290463863016e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (10.0, 2.0),
        9.95073771488337154e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (20.0, 2.0),
        9.98754668053816452e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (100.0, 2.0),
        9.99950007498750219e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1000.0, 2.0),
        9.99999500000749945e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (10000.0, 2.0),
        9.999999950000000739e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.0, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1e-100, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.001, 2.0),
        4.99646446697795041e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.01, 2.0),
        4.96464554479100486e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.1, 2.0),
        4.64732719207070087e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.325, 2.0),
        3.88014227253126233e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1.0, 2.0),
        2.11324865405187118e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1.5, 2.0),
        1.36196562445500540e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (2.0, 2.0),
        9.17517095361369836e-2,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (10.0, 2.0),
        4.92622851166284542e-3,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (20.0, 2.0),
        1.24533194618354849e-3,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (100.0, 2.0),
        4.99925012497812894e-5,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1000.0, 2.0),
        4.99999250001249998e-7,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (10000.0, 2.0),
        4.99999992500000125e-9,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1e-100, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.001, 2.0),
        4.99646446697795041e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.01, 2.0),
        4.96464554479100486e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.1, 2.0),
        4.64732719207070087e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.325, 2.0),
        3.88014227253126233e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1.0, 2.0),
        2.11324865405187118e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1.5, 2.0),
        1.36196562445500540e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-2.0, 2.0),
        9.17517095361369836e-02,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-10.0, 2.0),
        4.92622851166284542e-03,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-20.0, 2.0),
        1.24533194618354849e-03,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-100.0, 2.0),
        4.99925012497812894e-05,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1000.0, 2.0),
        4.99999250001249998e-07,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-10000.0, 2.0),
        4.99999992500000125e-09,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1e-100, 2.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.001, 2.0),
        5.00353553302204959e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.01, 2.0),
        5.03535445520899514e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.1, 2.0),
        5.35267280792929913e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.325, 2.0),
        6.11985772746873767e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1.0, 2.0),
        7.88675134594812882e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1.5, 2.0),
        8.63803437554499460e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-2.0, 2.0),
        9.08248290463863016e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-10.0, 2.0),
        9.95073771488337155e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-20.0, 2.0),
        9.98754668053816452e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-100.0, 2.0),
        9.99950007498750219e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1000.0, 2.0),
        9.99999500000749999e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-10000.0, 2.0),
        9.99999995000000075e-1,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.0, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1e-100, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.001, 300.0),
        5.00398609900942949e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.01, 300.0),
        5.03986033020559088e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.1, 300.0),
        5.39794441177768194e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (0.325, 300.0),
        6.27296201542523812e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1.0, 300.0),
        8.40941797784686861e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1.5, 300.0),
        9.32666983425369137e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (2.0, 300.0),
        9.76799239508425455e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (10.0, 300.0),
        1.00000000000000000e+00,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (20.0, 300.0),
        1.00000000000000000e+00,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (100.0, 300.0),
        1.00000000000000000e+00,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (1000.0, 300.0),
        1.00000000000000000e+00,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (10000.0, 300.0),
        1.00000000000000000e+00,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.0, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1e-100, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.001, 300.0),
        4.99601390099057051e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.01, 300.0),
        4.96013966979440912e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.1, 300.0),
        4.60205558822231806e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (0.325, 300.0),
        3.72703798457476188e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1.0, 300.0),
        1.59058202215313138e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1.5, 300.0),
        6.73330165746308628e-2,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (2.0, 300.0),
        2.32007604915745452e-2,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (10.0, 300.0),
        8.279313677e-21,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (20.0, 300.0),
        1.93159812815803978e-57,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (100.0, 300.0),
        1.02557519997736154e-232,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (1000.0, 300.0),
        0.00000000000000000e+00,
        0.0
    );
    check!(
        student_t_survival,
        (10000.0, 300.0),
        0.00000000000000000e+00,
        0.0
    );
    check!(
        student_t_cdf,
        (-1e-100, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.001, 300.0),
        4.99601390099057051e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.01, 300.0),
        4.96013966979440912e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.1, 300.0),
        4.60205558822231806e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-0.325, 300.0),
        3.72703798457476188e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1.0, 300.0),
        1.59058202215313138e-01,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-1.5, 300.0),
        6.73330165746308628e-02,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-2.0, 300.0),
        2.32007604915745452e-02,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-10.0, 300.0),
        8.279313675556272534e-21,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-20.0, 300.0),
        1.93159812815803978e-57,
        TEST_TOL6
    );
    check!(
        student_t_cdf,
        (-100.0, 300.0),
        1.02557519997736154e-232,
        TEST_TOL6
    );
    check!(student_t_cdf, (-1000.0, 300.0), 0.0, 0.0);
    check!(student_t_cdf, (-10000.0, 300.0), 0.0, 0.0);
    check!(
        student_t_survival,
        (-1e-100, 300.0),
        5.00000000000000000e-01,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.001, 300.0),
        5.00398609900942949e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.01, 300.0),
        5.03986033020559088e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.1, 300.0),
        5.39794441177768194e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-0.325, 300.0),
        6.27296201542523812e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1.0, 300.0),
        8.40941797784686862e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-1.5, 300.0),
        9.32666983425369137e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-2.0, 300.0),
        9.76799239508425455e-1,
        TEST_TOL6
    );
    check!(
        student_t_survival,
        (-10.0, 300.0),
        1.000000000000000000e0,
        TEST_TOL6
    );
    check!(student_t_survival, (-20.0, 300.0), 1.0, TEST_TOL6);
    check!(student_t_survival, (-100.0, 300.0), 1.0, TEST_TOL6);
    check!(student_t_survival, (-1000.0, 300.0), 1.0, TEST_TOL6);
    check!(student_t_survival, (-10000.0, 300.0), 1.0, TEST_TOL6);
}

/// Checks fdist.
#[test]
fn fdist() {
    check!(fisher_f_cdf, (0.0, 1.2, 1.3), 0.0, 0.0);
    check!(
        fisher_f_cdf,
        (1e-100, 1.2, 1.3),
        6.98194275525039002e-61,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.001, 1.2, 1.3),
        1.10608485860238564e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.01, 1.2, 1.3),
        4.38636757068313850e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.1, 1.2, 1.3),
        1.68242392712840734e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.325, 1.2, 1.3),
        3.14130045246195449e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.0, 1.2, 1.3),
        5.09630779074755253e-01,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.5, 1.2, 1.3),
        5.83998640641553852e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (2.0, 1.2, 1.3),
        6.34733581351938787e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (10.0, 1.2, 1.3),
        8.48446237879200975e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (20.0, 1.2, 1.3),
        9.00987726336875039e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (100.0, 1.2, 1.3),
        9.64489127047688435e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1000.0, 1.2, 1.3),
        9.92012051694116388e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (10000.0, 1.2, 1.3),
        9.98210862808842585e-1,
        TEST_TOL6
    );
    check!(fisher_f_survival, (0.0, 1.2, 1.3), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (1e-100, 1.2, 1.3), 1.0, TEST_TOL6);
    check!(
        fisher_f_survival,
        (0.001, 1.2, 1.3),
        9.88939151413976144e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.01, 1.2, 1.3),
        9.56136324293168615e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.1, 1.2, 1.3),
        8.31757607287159265e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.325, 1.2, 1.3),
        6.85869954753804551e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.0, 1.2, 1.3),
        4.90369220925244747e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.5, 1.2, 1.3),
        4.16001359358446148e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (2.0, 1.2, 1.3),
        3.65266418648061213e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10.0, 1.2, 1.3),
        1.51553762120799025e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (20.0, 1.2, 1.3),
        9.90122736631249612e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (100.0, 1.2, 1.3),
        3.55108729523115643e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1000.0, 1.2, 1.3),
        7.98794830588361109e-3,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10000.0, 1.2, 1.3),
        1.7891371911574145e-3,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (3.479082213465832574, 1.0, 4040712.0),
        0.93785072763723411967,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (3.002774644786533109, 1.0, 4040712.0),
        0.91687787379476055771,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (3.000854441173130827, 1.0, 4040712.0),
        0.91677930719813578619,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (3.000064021622133037, 1.0, 4040712.0),
        0.9167386972447996480399,
        TEST_TOL6
    );
    check!(fisher_f_cdf, (0.0, 500.0, 1.3), 0.0, 0.0);
    check!(fisher_f_cdf, (1e-100, 500.0, 1.3), 0.0, 0.0);
    check!(
        fisher_f_cdf,
        (0.001, 500.0, 1.3),
        9.83434460393304765e-141,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.01, 500.0, 1.3),
        1.45915624888550014e-26,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.1, 500.0, 1.3),
        5.89976509619688165e-4,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.325, 500.0, 1.3),
        6.86110486051542533e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.0, 500.0, 1.3),
        3.38475053806404615e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.5, 500.0, 1.3),
        4.52016245247457422e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (2.0, 500.0, 1.3),
        5.27339068937388798e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (10.0, 500.0, 1.3),
        8.16839628578413905e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (20.0, 500.0, 1.3),
        8.81784623056911406e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (100.0, 500.0, 1.3),
        9.58045057204221295e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1000.0, 500.0, 1.3),
        9.90585749380655275e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (10000.0, 500.0, 1.3),
        9.97891924831461387e-1,
        TEST_TOL6
    );
    check!(fisher_f_survival, (0.0, 500.0, 1.3), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (1e-100, 500.0, 1.3), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.001, 500.0, 1.3), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.01, 500.0, 1.3), 1.0, TEST_TOL6);
    check!(
        fisher_f_survival,
        (0.1, 500.0, 1.3),
        9.99410023490380312e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.325, 500.0, 1.3),
        9.31388951394845747e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.0, 500.0, 1.3),
        6.61524946193595385e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.5, 500.0, 1.3),
        5.47983754752542572e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (2.0, 500.0, 1.3),
        4.72660931062611202e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10.0, 500.0, 1.3),
        1.83160371421586096e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (20.0, 500.0, 1.3),
        1.18215376943088595e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (100.0, 500.0, 1.3),
        4.19549427957787016e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1000.0, 500.0, 1.3),
        9.41425061934473424e-3,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10000.0, 500.0, 1.3),
        2.10807516853862603e-3,
        TEST_TOL6
    );
    check!(fisher_f_cdf, (0.0, 1.2, 500.0), 0.0, 0.0);
    check!(
        fisher_f_cdf,
        (1e-100, 1.2, 500.0),
        8.23342055585482999e-61,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.001, 1.2, 500.0),
        1.30461496441289529e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.01, 1.2, 500.0),
        5.18324224608033294e-2,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.1, 1.2, 500.0),
        2.02235101716076289e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.325, 1.2, 500.0),
        3.90502983219393749e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.0, 1.2, 500.0),
        6.67656191574653619e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.5, 1.2, 500.0),
        7.75539230271467054e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (2.0, 1.2, 500.0),
        8.45209114904613705e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (10.0, 1.2, 500.0),
        9.99168017659120988e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (20.0, 1.2, 500.0),
        9.99998005738371669e-1,
        TEST_TOL6
    );
    check!(fisher_f_cdf, (100.0, 1.2, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (1000.0, 1.2, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (10000.0, 1.2, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.0, 1.2, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (1e-100, 1.2, 500.0), 1.0, TEST_TOL6);
    check!(
        fisher_f_survival,
        (0.001, 1.2, 500.0),
        9.86953850355871047e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.01, 1.2, 500.0),
        9.48167577539196671e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.1, 1.2, 500.0),
        7.97764898283923711e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (0.325, 1.2, 500.0),
        6.09497016780606251e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.0, 1.2, 500.0),
        3.32343808425346381e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.5, 1.2, 500.0),
        2.24460769728532946e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (2.0, 1.2, 500.0),
        1.54790885095386295e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10.0, 1.2, 500.0),
        8.3198234087901168e-4,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (20.0, 1.2, 500.0),
        1.99426162833131e-6,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (100.0, 1.2, 500.0),
        6.23302662288217117e-25,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1000.0, 1.2, 500.0),
        1.14328577259666930e-134,
        TEST_TOL6
    );
    check!(fisher_f_survival, (10000.0, 1.2, 500.0), 0.0, 0.0);
    check!(fisher_f_cdf, (0.0, 200.0, 500.0), 0.0, 0.0);
    check!(fisher_f_cdf, (1e-100, 200.0, 500.0), 0.0, 0.0);
    check!(
        fisher_f_cdf,
        (0.001, 200.0, 500.0),
        4.09325080403669893e-251,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.01, 200.0, 500.0),
        1.17894325419628688e-151,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.1, 200.0, 500.0),
        5.92430940796861258e-57,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (0.325, 200.0, 500.0),
        3.18220452357263554e-18,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.0, 200.0, 500.0),
        5.06746326121168266e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (1.5, 200.0, 500.0),
        9.99794175718712438e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_cdf,
        (2.0, 200.0, 500.0),
        9.99999999528236152e-1,
        TEST_TOL6
    );
    check!(fisher_f_cdf, (10.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (20.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (100.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (1000.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_cdf, (10000.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.0, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (1e-100, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.001, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.01, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(fisher_f_survival, (0.1, 200.0, 500.0), 1.0, TEST_TOL6);
    check!(
        fisher_f_survival,
        (0.325, 200.0, 500.0),
        9.99999999999999997e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.0, 200.0, 500.0),
        4.93253673878831734e-1,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (1.5, 200.0, 500.0),
        2.05824281287561795e-4,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (2.0, 200.0, 500.0),
        4.71763848371410786e-10,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (10.0, 200.0, 500.0),
        5.98048337181948436e-96,
        TEST_TOL6
    );
    check!(
        fisher_f_survival,
        (20.0, 200.0, 500.0),
        2.92099265879979502e-155,
        TEST_TOL6
    );
    check!(fisher_f_survival, (1000.0, 200.0, 500.0), 0.0, 0.0);
    check!(fisher_f_survival, (10000.0, 200.0, 500.0), 0.0, 0.0);
}

/// Checks finv.
#[test]
fn finv() {
    check!(fisher_f_quantile, (0.0, 1.2, 1.3), 0.0, 0.0);
    check!(
        fisher_f_quantile,
        (6.98194275525039002e-61, 1.2, 1.3),
        1e-100,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.10608485860238564e-2, 1.2, 1.3),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (4.38636757068313850e-2, 1.2, 1.3),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.68242392712840734e-1, 1.2, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (3.14130045246195449e-1, 1.2, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.09630779074755253e-01, 1.2, 1.3),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.83998640641553852e-1, 1.2, 1.3),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (6.34733581351938787e-1, 1.2, 1.3),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (8.48446237879200975e-1, 1.2, 1.3),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.00987726336875039e-1, 1.2, 1.3),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.64489127047688435e-1, 1.2, 1.3),
        100.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.92012051694116388e-1, 1.2, 1.3),
        1000.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.98210862808842585e-1, 1.2, 1.3),
        10000.0,
        TEST_TOL6
    );
    check!(fisher_f_inverse_survival, (1.0, 1.2, 1.3), 0.0, TEST_TOL6);
    check!(
        fisher_f_inverse_survival,
        (9.88939151413976144e-1, 1.2, 1.3),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (9.56136324293168615e-1, 1.2, 1.3),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (8.31757607287159265e-1, 1.2, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (6.85869954753804551e-1, 1.2, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.90369220925244747e-1, 1.2, 1.3),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.16001359358446148e-1, 1.2, 1.3),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (3.65266418648061213e-1, 1.2, 1.3),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.51553762120799025e-1, 1.2, 1.3),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (9.90122736631249612e-2, 1.2, 1.3),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (3.55108729523115643e-2, 1.2, 1.3),
        100.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (7.98794830588361109e-3, 1.2, 1.3),
        1000.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.7891371911574145e-3, 1.2, 1.3),
        10000.0,
        TEST_TOL6
    );
    check!(fisher_f_quantile, (0.0, 500.0, 1.3), 0.0, 0.0);
    check!(
        fisher_f_quantile,
        (9.83434460393304765e-141, 500.0, 1.3),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.45915624888550014e-26, 500.0, 1.3),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.89976509619688165e-4, 500.0, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (6.86110486051542533e-2, 500.0, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (3.38475053806404615e-1, 500.0, 1.3),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (4.52016245247457422e-1, 500.0, 1.3),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.27339068937388798e-1, 500.0, 1.3),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (8.16839628578413905e-1, 500.0, 1.3),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (8.81784623056911406e-1, 500.0, 1.3),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.58045057204221295e-1, 500.0, 1.3),
        100.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.90585749380655275e-1, 500.0, 1.3),
        1000.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.97891924831461387e-1, 500.0, 1.3),
        10000.0,
        TEST_TOL6
    );
    check!(fisher_f_inverse_survival, (1.0, 500.0, 1.3), 0.0, TEST_TOL6);
    check!(
        fisher_f_inverse_survival,
        (9.99410023490380312e-1, 500.0, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (9.31388951394845747e-1, 500.0, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (6.61524946193595385e-1, 500.0, 1.3),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (5.47983754752542572e-1, 500.0, 1.3),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.72660931062611202e-1, 500.0, 1.3),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.83160371421586096e-1, 500.0, 1.3),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.18215376943088595e-1, 500.0, 1.3),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.19549427957787016e-2, 500.0, 1.3),
        100.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (9.41425061934473424e-3, 500.0, 1.3),
        1000.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (2.10807516853862603e-3, 500.0, 1.3),
        10000.0,
        TEST_TOL6
    );
    check!(fisher_f_quantile, (0.0, 1.2, 500.0), 0.0, 0.0);
    check!(
        fisher_f_quantile,
        (8.23342055585482999e-61, 1.2, 500.0),
        1e-100,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.30461496441289529e-2, 1.2, 500.0),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.18324224608033294e-2, 1.2, 500.0),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (2.02235101716076289e-1, 1.2, 500.0),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (3.90502983219393749e-1, 1.2, 500.0),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (6.67656191574653619e-1, 1.2, 500.0),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (7.75539230271467054e-1, 1.2, 500.0),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (8.45209114904613705e-1, 1.2, 500.0),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.99168017659120988e-1, 1.2, 500.0),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.99998005738371669e-1, 1.2, 500.0),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.0, 1.2, 500.0),
        f64::INFINITY,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.0, 1.2, 500.0),
        f64::INFINITY,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.0, 1.2, 500.0),
        f64::INFINITY,
        TEST_TOL6
    );
    check!(fisher_f_inverse_survival, (1.0, 1.2, 500.0), 0.0, TEST_TOL6);
    check!(
        fisher_f_inverse_survival,
        (9.86953850355871047e-1, 1.2, 500.0),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (9.48167577539196671e-1, 1.2, 500.0),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (7.97764898283923711e-1, 1.2, 500.0),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (6.09497016780606251e-1, 1.2, 500.0),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (3.32343808425346381e-1, 1.2, 500.0),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (2.24460769728532946e-1, 1.2, 500.0),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.54790885095386295e-1, 1.2, 500.0),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (8.3198234087901168e-4, 1.2, 500.0),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.99426162833131e-6, 1.2, 500.0),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (6.23302662288217117e-25, 1.2, 500.0),
        100.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.14328577259666930e-134, 1.2, 500.0),
        1000.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (0.0, 1.2, 500.0),
        f64::INFINITY,
        0.0
    );
    check!(fisher_f_quantile, (0.0, 200.0, 500.0), 0.0, 0.0);
    check!(
        fisher_f_quantile,
        (4.09325080403669893e-251, 200.0, 500.0),
        0.001,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.17894325419628688e-151, 200.0, 500.0),
        0.01,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.92430940796861258e-57, 200.0, 500.0),
        0.1,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (3.18220452357263554e-18, 200.0, 500.0),
        0.325,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (5.06746326121168266e-1, 200.0, 500.0),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (9.99794175718712438e-1, 200.0, 500.0),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_quantile,
        (1.0, 200.0, 500.0),
        f64::INFINITY,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (1.0, 200.0, 500.0),
        0.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.93253673878831734e-1, 200.0, 500.0),
        1.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (2.05824281287561795e-4, 200.0, 500.0),
        1.5,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (4.71763848371410786e-10, 200.0, 500.0),
        2.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (5.98048337181948436e-96, 200.0, 500.0),
        10.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (2.92099265879979502e-155, 200.0, 500.0),
        20.0,
        TEST_TOL6
    );
    check!(
        fisher_f_inverse_survival,
        (0.0, 200.0, 500.0),
        f64::INFINITY,
        0.0
    );
    check!(
        fisher_f_quantile,
        (0.95, 1.0, 261.0),
        3.8773340322508720313e+00,
        TEST_TOL3
    );
}

/// Checks gamma.
#[test]
fn gamma() {
    check!(gamma_cdf, (0.0, 1.0, 1.0), 0.0, 0.0);
    check!(gamma_cdf, (1e-100, 1.0, 1.0), 1e-100, TEST_TOL6);
    check!(
        gamma_cdf,
        (0.001, 1.0, 1.0),
        9.99500166625008332e-4,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.01, 1.0, 1.0),
        9.95016625083194643e-3,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.1, 1.0, 1.0),
        9.51625819640404268e-2,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.325, 1.0, 1.0),
        2.77472646357927811e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.0, 1.0, 1.0),
        6.32120558828557678e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.5, 1.0, 1.0),
        7.76869839851570171e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (2.0, 1.0, 1.0),
        8.64664716763387308e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (10.0, 1.0, 1.0),
        9.99954600070237515e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (20.0, 1.0, 1.0),
        9.99999997938846378e-1,
        TEST_TOL6
    );
    check!(gamma_cdf, (100.0, 1.0, 1.0), 1e0, TEST_TOL6);
    check!(gamma_cdf, (1000.0, 1.0, 1.0), 1e0, TEST_TOL6);
    check!(gamma_cdf, (10000.0, 1.0, 1.0), 1e0, TEST_TOL6);
    check!(gamma_survival, (0.0, 1.0, 1.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (1e-100, 1.0, 1.0), 1.0, TEST_TOL6);
    check!(
        gamma_survival,
        (0.001, 1.0, 1.0),
        9.99000499833374992e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.01, 1.0, 1.0),
        9.90049833749168054e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.1, 1.0, 1.0),
        9.04837418035959573e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.325, 1.0, 1.0),
        7.22527353642072189e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1.0, 1.0, 1.0),
        3.67879441171442322e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1.5, 1.0, 1.0),
        2.23130160148429829e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (2.0, 1.0, 1.0),
        1.35335283236612692e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (10.0, 1.0, 1.0),
        4.53999297624848515e-5,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (20.0, 1.0, 1.0),
        2.06115362243855783e-9,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (100.0, 1.0, 1.0),
        3.72007597602083596e-44,
        TEST_TOL6
    );
    check!(gamma_survival, (1000.0, 1.0, 1.0), 0.0, 0.0);
    check!(gamma_survival, (10000.0, 1.0, 1.0), 0.0, 0.0);
    check!(gamma_cdf, (0.0, 1.0, 10.0), 0.0, 0.0);
    check!(gamma_cdf, (1e-100, 1.0, 10.0), 1e-101, TEST_TOL6);
    check!(
        gamma_cdf,
        (0.001, 1.0, 10.0),
        9.99950001666625001e-5,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.01, 1.0, 10.0),
        9.99500166625008332e-4,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.1, 1.0, 10.0),
        9.95016625083194643e-3,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.325, 1.0, 10.0),
        3.19775501686939529e-2,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.0, 1.0, 10.0),
        9.51625819640404268e-2,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.5, 1.0, 10.0),
        1.39292023574942193e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (2.0, 1.0, 10.0),
        1.81269246922018141e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (10.0, 1.0, 10.0),
        6.32120558828557678e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (20.0, 1.0, 10.0),
        8.64664716763387308e-1,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (100.0, 1.0, 10.0),
        9.99954600070237515e-1,
        TEST_TOL6
    );
    check!(gamma_cdf, (1000.0, 1.0, 10.0), 1e0, TEST_TOL6);
    check!(gamma_cdf, (10000.0, 1.0, 10.0), 1e0, TEST_TOL6);
    check!(gamma_survival, (0.0, 1.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (1e-100, 1.0, 10.0), 1.0, TEST_TOL6);
    check!(
        gamma_survival,
        (0.001, 1.0, 10.0),
        9.99900004999833337e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.01, 1.0, 10.0),
        9.99000499833374992e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.1, 1.0, 10.0),
        9.90049833749168054e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (0.325, 1.0, 10.0),
        9.68022449831306047e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1.0, 1.0, 10.0),
        9.04837418035959573e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1.5, 1.0, 10.0),
        8.60707976425057807e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (2.0, 1.0, 10.0),
        8.18730753077981859e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (10.0, 1.0, 10.0),
        3.67879441171442322e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (20.0, 1.0, 10.0),
        1.35335283236612692e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (100.0, 1.0, 10.0),
        4.53999297624848515e-5,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1000.0, 1.0, 10.0),
        3.72007597602083596e-44,
        TEST_TOL6
    );
    check!(gamma_survival, (10000.0, 1.0, 10.0), 0.0, 0.0);
    check!(gamma_cdf, (0.0, 17.0, 10.0), 0e0, 0.0);
    check!(gamma_cdf, (1e-100, 17.0, 10.0), 0e0, 0.0);
    check!(
        gamma_cdf,
        (0.001, 17.0, 10.0),
        2.81119174040422844e-83,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.01, 17.0, 10.0),
        2.80880324651985887e-66,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.1, 17.0, 10.0),
        2.78502998087492130e-49,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (0.325, 17.0, 10.0),
        1.37283653245125844e-40,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.0, 17.0, 10.0),
        2.55811932292544243e-32,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1.5, 17.0, 10.0),
        2.40420441175422372e-29,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (2.0, 17.0, 10.0),
        3.05092926217898577e-27,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (10.0, 17.0, 10.0),
        1.094920130378183e-15,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (20.0, 17.0, 10.0),
        5.60605096173161688e-11,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (100.0, 17.0, 10.0),
        2.70416097848011280e-2,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (1000.0, 17.0, 10.0),
        1.000000000000000000e0,
        TEST_TOL6
    );
    check!(
        gamma_cdf,
        (10000.0, 17.0, 10.0),
        1.000000000000000000e0,
        TEST_TOL6
    );
    check!(gamma_survival, (0.0, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (1e-100, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (0.001, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (0.01, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (0.1, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (0.325, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (1.0, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (1.5, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(gamma_survival, (2.0, 17.0, 10.0), 1.0, TEST_TOL6);
    check!(
        gamma_survival,
        (10.0, 17.0, 10.0),
        9.99999999999998905e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (20.0, 17.0, 10.0),
        9.99999999943939490e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (100.0, 17.0, 10.0),
        9.72958390215198872e-1,
        TEST_TOL6
    );
    check!(
        gamma_survival,
        (1000.0, 17.0, 10.0),
        2.11200951633948570e-25,
        TEST_TOL6
    );
    check!(gamma_survival, (10000.0, 17.0, 10.0), 0.0, 0.0);
}

/// Checks chisq.
#[test]
fn chisq() {
    check!(chi_squared_cdf, (0.0, 13.0), 0.0, 0.0);
    check!(chi_squared_cdf, (1e-100, 13.0), 0.0, 0.0);
    check!(
        chi_squared_cdf,
        (0.001, 13.0),
        1.86631102655845996e-25,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (0.01, 13.0),
        5.87882248504529790e-19,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (0.1, 13.0),
        1.78796983358555410e-12,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (0.325, 13.0),
        3.44611313779905183e-9,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (1.0, 13.0),
        3.83473473513595154e-6,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (1.5, 13.0),
        4.31718389201041932e-5,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (2.0, 13.0),
        2.26250084656047180e-4,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (10.0, 13.0),
        3.06065632019251110e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (20.0, 13.0),
        9.04789743921908487e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_cdf,
        (100.0, 13.0),
        9.99999999999998341e-1,
        TEST_TOL6
    );
    check!(chi_squared_cdf, (1000.0, 13.0), 1e0, TEST_TOL6);
    check!(chi_squared_cdf, (10000.0, 13.0), 1e0, TEST_TOL6);
    check!(chi_squared_survival, (0.0, 13.0), 1e0, TEST_TOL6);
    check!(chi_squared_survival, (1e-100, 13.0), 1e0, TEST_TOL6);
    check!(chi_squared_survival, (0.001, 13.0), 1e0, TEST_TOL6);
    check!(
        chi_squared_survival,
        (0.01, 13.0),
        9.99999999999999999e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (0.1, 13.0),
        9.99999999998212030e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (0.325, 13.0),
        9.99999996553886862e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (1.0, 13.0),
        9.99996165265264864e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (1.5, 13.0),
        9.99956828161079896e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (2.0, 13.0),
        9.99773749915343953e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (10.0, 13.0),
        6.93934367980748890e-1,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (20.0, 13.0),
        9.52102560780915127e-2,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (100.0, 13.0),
        1.65902608070858809e-15,
        TEST_TOL6
    );
    check!(
        chi_squared_survival,
        (1000.0, 13.0),
        1.74851191544860225e-205,
        TEST_TOL6
    );
    check!(chi_squared_survival, (10000.0, 13.0), 0.0, 0.0);
}

/// Checks beta.
#[test]
fn beta() {
    check!(beta_cdf, (0.0, 1.2, 1.3), 0.0, 0.0);
    check!(
        beta_cdf,
        (1e-100, 1.2, 1.3),
        1.34434944656489596e-120,
        TEST_TOL6
    );
    check!(
        beta_cdf,
        (0.001, 1.2, 1.3),
        3.37630042504535813e-4,
        TEST_TOL6
    );
    check!(
        beta_cdf,
        (0.01, 1.2, 1.3),
        5.34317264038929473e-3,
        TEST_TOL6
    );
    check!(beta_cdf, (0.1, 1.2, 1.3), 8.33997828306748346e-2, TEST_TOL6);
    check!(
        beta_cdf,
        (0.325, 1.2, 1.3),
        3.28698654180583916e-1,
        TEST_TOL6
    );
    check!(beta_cdf, (0.5, 1.2, 1.3), 5.29781429451299081e-1, TEST_TOL6);
    check!(beta_cdf, (0.9, 1.2, 1.3), 9.38529397224430659e-1, TEST_TOL6);
    check!(
        beta_cdf,
        (0.99, 1.2, 1.3),
        9.96886438341254380e-1,
        TEST_TOL6
    );
    check!(
        beta_cdf,
        (0.999, 1.2, 1.3),
        9.99843792833067634e-1,
        TEST_TOL6
    );
    check!(beta_cdf, (1.0, 1.2, 1.3), 1.0, TEST_TOL6);
    check!(beta_survival, (0.0, 1.2, 1.3), 1.0, 0.0);
    check!(beta_survival, (1e-100, 1.2, 1.3), 1e0, TEST_TOL6);
    check!(
        beta_survival,
        (0.001, 1.2, 1.3),
        9.99662369957495464e-1,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.01, 1.2, 1.3),
        9.94656827359610705e-1,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.1, 1.2, 1.3),
        9.16600217169325165e-1,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.325, 1.2, 1.3),
        6.71301345819416084e-1,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.5, 1.2, 1.3),
        4.70218570548700919e-1,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.9, 1.2, 1.3),
        6.14706027755693408e-2,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.99, 1.2, 1.3),
        3.11356165874561958e-3,
        TEST_TOL6
    );
    check!(
        beta_survival,
        (0.999, 1.2, 1.3),
        1.56207166932365759e-4,
        TEST_TOL6
    );
    check!(beta_survival, (1.0, 1.2, 1.3), 0.0, TEST_TOL6);
}

/// Checks betainv.
#[test]
fn betainv() {
    check!(beta_quantile, (0.0, 1.2, 1.3), 0.0, 0.0);
    check!(
        beta_quantile,
        (1.34434944656489596e-120, 1.2, 1.3),
        1e-100,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (3.37630042504535813e-4, 1.2, 1.3),
        0.001,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (5.34317264038929473e-3, 1.2, 1.3),
        0.01,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (8.33997828306748346e-2, 1.2, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (3.28698654180583916e-1, 1.2, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (5.29781429451299081e-1, 1.2, 1.3),
        0.5,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (9.38529397224430659e-1, 1.2, 1.3),
        0.9,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (9.96886438341254380e-1, 1.2, 1.3),
        0.99,
        TEST_TOL6
    );
    check!(
        beta_quantile,
        (9.99843792833067634e-1, 1.2, 1.3),
        0.999,
        TEST_TOL6
    );
    check!(beta_quantile, (1.0, 1.2, 1.3), 1.0, TEST_TOL6);
    check!(beta_inverse_survival, (1.0, 1.2, 1.3), 0.0, 0.0);
    check!(beta_inverse_survival, (1e0, 1.2, 1.3), 0.0, TEST_TOL6);
    check!(
        beta_inverse_survival,
        (9.99662369957495464e-1, 1.2, 1.3),
        0.001,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (9.94656827359610705e-1, 1.2, 1.3),
        0.01,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (9.16600217169325165e-1, 1.2, 1.3),
        0.1,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (6.71301345819416084e-1, 1.2, 1.3),
        0.325,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (4.70218570548700919e-1, 1.2, 1.3),
        0.5,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (6.14706027755693408e-2, 1.2, 1.3),
        0.9,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (3.11356165874561958e-3, 1.2, 1.3),
        0.99,
        TEST_TOL6
    );
    check!(
        beta_inverse_survival,
        (1.56207166932365759e-4, 1.2, 1.3),
        0.999,
        TEST_TOL6
    );
    check!(beta_inverse_survival, (0.0, 1.2, 1.3), 1.0, TEST_TOL6);
    check!(
        beta_quantile,
        (0.025, 2133.0, 7868.0),
        0.20530562929915865457928654,
        TEST_TOL6
    );
}

/// Checks gammainv.
#[test]
fn gammainv() {
    check!(gamma_quantile, (0.0, 1.0, 1.0), 0.0, 0.0);
    check!(gamma_quantile, (1e-100, 1.0, 1.0), 1e-100, TEST_TOL6);
    check!(
        gamma_quantile,
        (9.99500166625008332e-4, 1.0, 1.0),
        0.001,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (9.95016625083194643e-3, 1.0, 1.0),
        0.01,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (9.51625819640404268e-2, 1.0, 1.0),
        0.1,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (2.77472646357927811e-1, 1.0, 1.0),
        0.325,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (6.32120558828557678e-1, 1.0, 1.0),
        1.0,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (7.76869839851570171e-1, 1.0, 1.0),
        1.5,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (8.64664716763387308e-1, 1.0, 1.0),
        2.0,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (9.99954600070237515e-1, 1.0, 1.0),
        10.0,
        TEST_TOL6
    );
    check!(
        gamma_quantile,
        (9.99999997938846378e-1, 1.0, 1.0),
        20.0,
        100.0 * TEST_TOL6
    );
    check!(gamma_quantile, (1.0, 1.0, 1.0), f64::INFINITY, 0.0);
    check!(
        gamma_quantile,
        (0.1, 11.887411491530846, 1.0),
        7.73788447848618e+00,
        TEST_TOL1
    );
    check!(gamma_inverse_survival, (0.0, 1.0, 1.0), f64::INFINITY, 0.0);
    check!(
        gamma_inverse_survival,
        (2.06115362243855783e-9, 1.0, 1.0),
        20.0,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (4.53999297624848515e-5, 1.0, 1.0),
        10.0,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (1.35335283236612692e-1, 1.0, 1.0),
        2.0,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (2.23130160148429829e-1, 1.0, 1.0),
        1.5,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (3.67879441171442322e-1, 1.0, 1.0),
        1.0,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (7.22527353642072189e-1, 1.0, 1.0),
        0.325,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (9.04837418035959573e-1, 1.0, 1.0),
        0.1,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (9.90049833749168054e-1, 1.0, 1.0),
        0.01,
        TEST_TOL6
    );
    check!(
        gamma_inverse_survival,
        (9.99000499833374992e-1, 1.0, 1.0),
        0.001,
        TEST_TOL6
    );
    check!(gamma_inverse_survival, (1.0, 1.0, 1.0), 0.0, 0.0);
}

/// Checks chisqinv.
#[test]
fn chisqinv() {
    check!(chi_squared_quantile, (0.0, 13.0), 0.0, 0.0);
    check!(
        chi_squared_quantile,
        (1.86631102655845996e-25, 13.0),
        0.001,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (5.87882248504529790e-19, 13.0),
        0.01,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (1.78796983358555410e-12, 13.0),
        0.1,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (3.44611313779905183e-9, 13.0),
        0.325,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (3.83473473513595154e-6, 13.0),
        1.0,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (4.31718389201041932e-5, 13.0),
        1.5,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (2.26250084656047180e-4, 13.0),
        2.0,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (3.06065632019251110e-1, 13.0),
        10.0,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (9.04789743921908487e-1, 13.0),
        20.0,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (9.99999999999998341e-1, 13.0),
        100.0,
        0.01
    );
    check!(chi_squared_quantile, (1e0, 13.0), f64::INFINITY, 0.0);
    check!(
        chi_squared_quantile,
        (1.93238145206123590e-01, 1.5),
        0.211980092931799521729407,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (4.83e-8, 19.19),
        1.632280186860266704532868343,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (0.05, 1263131.0),
        1260517.771133388726131469059,
        TEST_TOL6
    );
    check!(
        chi_squared_quantile,
        (0.05, 2526262.0),
        2522565.864973351096735720202,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (0.0, 13.0),
        f64::INFINITY,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (1.65902608070858809e-15, 13.0),
        100.0,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (9.52102560780915127e-2, 13.0),
        20.0,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (6.93934367980748892e-1, 13.0),
        10.0,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (9.99773749915343954e-1, 13.0),
        2.0,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (9.99956828161079894e-1, 13.0),
        1.5,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (9.99996165265264863e-1, 13.0),
        1.0,
        TEST_TOL6
    );
    check!(
        chi_squared_inverse_survival,
        (9.99999996553886862e-1, 13.0),
        0.325,
        1e-6
    );
    check!(
        chi_squared_inverse_survival,
        (9.99999999998212031e-1, 13.0),
        0.1,
        1e-5
    );
    check!(chi_squared_inverse_survival, (1.0, 13.0), 0.0, 0.0);
}

/// Checks tdistinv.
#[test]
fn tdistinv() {
    check!(student_t_quantile, (0.5, 1.0), 0.0, TEST_TOL6);
    check!(
        student_t_quantile,
        (5.00318309780080559e-1, 1.0),
        0.001,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (5.03182992764908255e-1, 1.0),
        0.01,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (5.31725517430553569e-1, 1.0),
        0.1,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (6.00023120032852123e-1, 1.0),
        0.325,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (0.75000000000000000e0, 1.0),
        1.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (8.12832958189001183e-1, 1.0),
        1.5,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (8.52416382349566726e-1, 1.0),
        2.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.68274482569446430e-1, 1.0),
        10.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.84097748743823625e-1, 1.0),
        20.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.96817007235091745e-1, 1.0),
        100.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.99681690219919441e-1, 1.0),
        1000.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.99968169011487724e-1, 1.0),
        10000.0,
        TEST_TOL6
    );
    check!(student_t_inverse_survival, (0.5, 1.0), 0.0, TEST_TOL6);
    check!(
        student_t_inverse_survival,
        (4.99681690219919441e-1, 1.0),
        0.001,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (4.96817007235091745e-1, 1.0),
        0.01,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (4.68274482569446430e-1, 1.0),
        0.1,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.99976879967147876e-1, 1.0),
        0.325,
        TEST_TOL6
    );
    check!(student_t_inverse_survival, (2.5e-1, 1.0), 1.0, TEST_TOL6);
    check!(
        student_t_inverse_survival,
        (1.87167041810998816e-1, 1.0),
        1.5,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.47583617650433274e-1, 1.0),
        2.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.17255174305535695e-2, 1.0),
        10.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.59022512561763752e-2, 1.0),
        20.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.18299276490825515e-3, 1.0),
        100.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.18309780080558939e-4, 1.0),
        1000.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.18309885122757724e-5, 1.0),
        10000.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.99681690219919441e-1, 1.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.96817007235091744e-1, 1.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.68274482569446430e-1, 1.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.99976879967147876e-1, 1.0),
        -0.325,
        TEST_TOL6
    );
    check!(student_t_quantile, (0.25, 1.0), -1.0, TEST_TOL6);
    check!(
        student_t_quantile,
        (1.87167041810998816e-1, 1.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.47583617650433274e-1, 1.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.17255174305535695e-2, 1.0),
        -10.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.59022512561763751e-2, 1.0),
        -20.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.18299276490825514e-3, 1.0),
        -100.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.18309780080558938e-4, 1.0),
        -1000.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.18309885122757724e-5, 1.0),
        -10000.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.00318309780080559e-1, 1.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.03182992764908255e-1, 1.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.31725517430553570e-1, 1.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (6.00023120032852124e-1, 1.0),
        -0.325,
        TEST_TOL6
    );
    check!(student_t_inverse_survival, (7.5e-1, 1.0), -1.0, TEST_TOL6);
    check!(
        student_t_inverse_survival,
        (8.12832958189001184e-1, 1.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (8.52416382349566726e-1, 1.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.68274482569446430e-1, 1.0),
        -10.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.84097748743823625e-1, 1.0),
        -20.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.96817007235091745e-1, 1.0),
        -100.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.99681690219919441e-1, 1.0),
        -1000.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.99968169011487724e-1, 1.0),
        -10000.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.99646446697795041e-01, 2.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.96464554479100486e-01, 2.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.64732719207070087e-01, 2.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.88014227253126233e-01, 2.0),
        -0.325,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (2.11324865405187118e-01, 2.0),
        -1.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.36196562445500540e-01, 2.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.17517095361369836e-02, 2.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.92622851166284542e-03, 2.0),
        -10.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.24533194618354849e-03, 2.0),
        -20.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.99925012497812894e-05, 2.0),
        -100.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.99999250001249998e-07, 2.0),
        -1000.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.99999992500000125e-09, 2.0),
        -10000.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.00353553302204959e-1, 2.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.03535445520899514e-1, 2.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.35267280792929913e-1, 2.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (6.11985772746873767e-1, 2.0),
        -0.325,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (7.88675134594812882e-1, 2.0),
        -1.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (8.63803437554499460e-1, 2.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.08248290463863016e-1, 2.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.95073771488337155e-1, 2.0),
        -10.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.98754668053816452e-1, 2.0),
        -20.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.99950007498750219e-1, 2.0),
        -100.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.99999500000749999e-1, 2.0),
        -1000.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.99999995000000075e-1, 2.0),
        -10000.0,
        1e-6
    );
    check!(
        student_t_quantile,
        (5.00000000000000000e-01, 300.0),
        0.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (5.00398609900942949e-01, 300.0),
        0.001,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (5.03986033020559088e-01, 300.0),
        0.01,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (5.39794441177768194e-01, 300.0),
        0.1,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (6.27296201542523812e-01, 300.0),
        0.325,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (8.40941797784686861e-01, 300.0),
        1.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.32666983425369137e-01, 300.0),
        1.5,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (9.76799239508425455e-01, 300.0),
        2.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.00000000000000000e+00, 300.0),
        f64::INFINITY,
        0.0
    );
    check!(
        student_t_inverse_survival,
        (5.00000000000000000e-01, 300.0),
        0.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (4.99601390099057051e-1, 300.0),
        0.001,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (4.96013966979440912e-1, 300.0),
        0.01,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (4.60205558822231806e-1, 300.0),
        0.1,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (3.72703798457476188e-1, 300.0),
        0.325,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.59058202215313138e-1, 300.0),
        1.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (6.73330165746308628e-2, 300.0),
        1.5,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (2.32007604915745452e-2, 300.0),
        2.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (8.279313677e-21, 300.0),
        10.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.93159812815803978e-57, 300.0),
        20.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.02557519997736154e-232, 300.0),
        100.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (0.00000000000000000e+00, 300.0),
        f64::INFINITY,
        0.0
    );
    check!(
        student_t_quantile,
        (4.99601390099057051e-01, 300.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.96013966979440912e-01, 300.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (4.60205558822231806e-01, 300.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (3.72703798457476188e-01, 300.0),
        -0.325,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.59058202215313138e-01, 300.0),
        -1.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (6.73330165746308628e-02, 300.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (2.32007604915745452e-02, 300.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (8.279313675556272534e-21, 300.0),
        -10.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.93159812815803978e-57, 300.0),
        -20.0,
        TEST_TOL6
    );
    check!(
        student_t_quantile,
        (1.02557519997736154e-232, 300.0),
        -100.0,
        TEST_TOL6
    );
    check!(student_t_quantile, (0.0, 300.0), f64::NEG_INFINITY, 0.0);
    check!(
        student_t_inverse_survival,
        (5.00398609900942949e-1, 300.0),
        -0.001,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.03986033020559088e-1, 300.0),
        -0.01,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (5.39794441177768194e-1, 300.0),
        -0.1,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (6.27296201542523812e-1, 300.0),
        -0.325,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (8.40941797784686862e-1, 300.0),
        -1.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.32666983425369137e-1, 300.0),
        -1.5,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (9.76799239508425455e-1, 300.0),
        -2.0,
        TEST_TOL6
    );
    check!(
        student_t_inverse_survival,
        (1.000000000000000000e0, 300.0),
        f64::NEG_INFINITY,
        TEST_TOL6
    );
}
