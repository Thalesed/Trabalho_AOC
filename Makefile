CC ?= cc
CPPFLAGS ?= -I.
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror
LDLIBS ?= -lm
USE_FLOAT128 ?= 1
ifeq ($(USE_FLOAT128),1)
CPPFLAGS += -DTAKUM_USE_FLOAT128
LDLIBS += -lquadmath
endif
BUILD := build
CORE := Takum.c src/rvv_takum.c
HEADERS := takum.h src/rvv_takum.h src/tk_wide.h
LIBTAKUM ?= ../reference-libtakum
.PHONY: all demo test sanitize oracle precision clean
all: $(BUILD)/demo_scalar $(BUILD)/demo_rvv
$(BUILD):
	mkdir -p $@
$(BUILD)/demo_scalar: Takum.c takum.h src/tk_wide.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) Takum.c $(LDLIBS) -o $@
$(BUILD)/demo_rvv: $(CORE) $(HEADERS) examples/demo_rvv.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DTAKUM_NO_MAIN $(CORE) examples/demo_rvv.c $(LDLIBS) -o $@
$(BUILD)/test_takum: $(CORE) $(HEADERS) tests/test_takum.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DTAKUM_NO_MAIN $(CORE) tests/test_takum.c $(LDLIBS) -o $@
demo: all
	./$(BUILD)/demo_scalar
	./$(BUILD)/demo_rvv 8
	./$(BUILD)/demo_rvv 16
	./$(BUILD)/demo_rvv 32
test: $(BUILD)/test_takum
	./$(BUILD)/test_takum
$(BUILD)/test_float128: Takum.c $(HEADERS) tests/test_float128.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DTAKUM_NO_MAIN Takum.c tests/test_float128.c $(LDLIBS) -lquadmath -o $@
precision: $(BUILD)/test_float128
	./$(BUILD)/test_float128
sanitize: | $(BUILD)
	$(CC) $(CPPFLAGS) -std=c11 -O1 -g -Wall -Wextra -Werror -fsanitize=address,undefined -fno-omit-frame-pointer -DTAKUM_NO_MAIN $(CORE) tests/test_takum.c $(LDLIBS) -o $(BUILD)/test_sanitize
	ASAN_OPTIONS=detect_leaks=0 ./$(BUILD)/test_sanitize
oracle: | $(BUILD)
	$(MAKE) -C $(LIBTAKUM) libtakum.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -DLIBTAKUM_HEADER='"$(abspath $(LIBTAKUM))/takum.h"' -DTAKUM_NO_MAIN Takum.c tests/test_oracle.c $(LIBTAKUM)/libtakum.a $(LDLIBS) -o $(BUILD)/test_oracle
	./$(BUILD)/test_oracle
clean:
	rm -rf $(BUILD)
