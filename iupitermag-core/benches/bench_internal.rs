use criterion::{Criterion, criterion_group, criterion_main};
use iupitermag::{field::Field, internal::InternalField};
use std::{f64::consts::PI, hint::black_box};

fn criterion_benchmark(c: &mut Criterion) {
    let internal_field = InternalField::new("JRM33", None, None, None);
    c.bench_function("Internal Field", |b| {
        b.iter(|| {
            _ = internal_field.calc_field(
                black_box(10.),
                black_box(PI * 0.25),
                black_box(PI * 0.3),
            );
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
