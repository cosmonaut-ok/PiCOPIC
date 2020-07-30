#include <gtest/gtest.h>
#include "msg.hpp"

#define LOGURU_WITH_STREAMS 1

TEST(Message, MSG)
{
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  MSG("Message");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(output, "Message\n");
  ASSERT_EQ(error, "");
}

TEST(Message, MSG_FIXME)
{
  // unfixme
#undef FIXME
#define FIXME false
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  MSG_FIXME("Fixme");
  std::string output_u = testing::internal::GetCapturedStdout();
  std::string error_u = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error_u, "");
  ASSERT_EQ(output_u, "");

  // fixme
#undef FIXME
#define FIXME true
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  MSG_FIXME("Fixme");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error, "FIXME: Fixme.\n");
  ASSERT_EQ(output, "");
}
