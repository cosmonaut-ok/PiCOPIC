#include <gtest/gtest.h>
#include "msg.hpp"


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

TEST(Message, LOG_CRIT)
{

  ASSERT_DEATH( {
      LOG_CRIT("Critical", 1);
    }, "ERROR: Critical. Can not continue. Exiting.");
}

TEST(Message, LOG_ERR)
{
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  LOG_ERR("Error");
  std::string error = testing::internal::GetCapturedStderr();
  std::string output = testing::internal::GetCapturedStdout();
  ASSERT_EQ(error, "ERROR: Error.\n");
  ASSERT_EQ(output, "");
}

TEST(Message, LOG_WARN)
{
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  LOG_WARN("Warning");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error, "WARNING: Warning.\n");
  ASSERT_EQ(output, "");
}

TEST(Message, LOG_INFO)
{
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  LOG_INFO("Info");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error, "INFO: Info.\n");
  ASSERT_EQ(output, "");
}

TEST(Message, LOG_DBG)
{
  // undebug mode
#define DEBUG false
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  LOG_DBG("Debug");
  std::string output_u = testing::internal::GetCapturedStdout();
  std::string error_u = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error_u, "");
  ASSERT_EQ(output_u, "");
  
  // debug mode
#define DEBUG true
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  LOG_DBG("Debug");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error, "DEBUG: Debug.\n");
  ASSERT_EQ(output, "");

}

TEST(Message, MSG_FIXME)
{
  // unfixme
#define FIXME false
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  MSG_FIXME("Fixme");
  std::string output_u = testing::internal::GetCapturedStdout();
  std::string error_u = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error_u, "");
  ASSERT_EQ(output_u, "");

  // fixme
#define FIXME true
  testing::internal::CaptureStderr();
  testing::internal::CaptureStdout();
  MSG_FIXME("Fixme");
  std::string output = testing::internal::GetCapturedStdout();
  std::string error = testing::internal::GetCapturedStderr();
  ASSERT_EQ(error, "FIXME: Fixme.\n");
  ASSERT_EQ(output, "");
}
